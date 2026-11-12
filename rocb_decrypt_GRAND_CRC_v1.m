function [c_decoded, M_out, valid, abandoned, n_guesses, info] = ...
    rocb_decrypt_GRAND_CRC_v1(key_hex, N, C_in, T_in, taglen, n, k, max_query, CRC)
%ROCB_DECRYPT_GRAND_CRC_V1  Joint CRC-R-OCB decryption/decoding using GRAND.
%
% This implementation follows the ROCB decryption pseudocode:
%   1. Split the noisy input into ciphertext blocks and tag.
%   2. Compute L = E_K(0^n), H = E_K(N xor L), and Z[i].
%   3. Decode each ciphertext block by GRAND membership tests on padding.
%   4. Use the special last-block GRAND rule with L*x^{-1}.
%   5. Decode the noisy tag with GRAND when tau = n, then verify it.
%
% Inputs
%   key_hex   : AES key as hex string, used by Cipher/InvCipher.
%   N         : nonce, either n-bit vector or n/4-character hex string.
%   C_in      : noisy ciphertext bits, without the tag.
%   T_in      : noisy tag bits. The pseudocode supports GRAND tag decoding
%               only when taglen = n.
%   taglen    : tag length in bits.
%   n         : block length in bits, e.g., 128.
%   k         : information length in bits, e.g., 112, 116, 120.
%   max_query : maximum GRAND guesses per block.
%   P         : fixed redundancy padding of length n-k bits. Optional;
%               default is all zeros.
%
% Outputs
%   c_decoded : decoded/corrected ciphertext bits, without the tag.
%   M_out     : decoded message bits, or [] if invalid/abandoned.
%   valid     : true if all blocks decoded and the tag verified.
%   abandoned : per-ciphertext-block abandonment indicators.
%   n_guesses : per-ciphertext-block GRAND query counts.
%   info      : debug information, including tag GRAND statistics.
%
%
% Regular data-block membership:
%       H_crc * D' = 0  over GF(2)
%
% Special data-block membership:
%       complement bit k+1 of D, then test CRC syndrome.
%
% The special-block marker is therefore the CRC analogue of P_tilde.
% The tag decoding rule remains the same as in the zero-padding R-OCB
% implementation because the R-OCB authentication tag itself is unchanged.

    % ---------------------- input checks ----------------------
    validateattributes(n, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(k, {'numeric'}, {'scalar','integer','positive','<',n});
    validateattributes(taglen, {'numeric'}, {'scalar','integer','positive','<=',n});

    if mod(n,8) ~= 0
        error('n must be a multiple of 8 because AES/Cipher is byte based.');
    end

    if ~isstruct(CRC) || ~isfield(CRC,'H')
        error('CRC must be a struct containing CRC.H.');
    end

    H_crc = double(CRC.H);

    % make_CRC_GH used in the GRAND codebase returns H as (n-k)-by-n.
    if size(H_crc,2) ~= n
        error('CRC.H must have n=%d columns.', n);
    end

    if isempty(C_in) || mod(length(C_in),n) ~= 0
        error('C_in length must be a positive multiple of n.');
    end

    if length(T_in) ~= taglen
        error('T_in must contain exactly taglen bits.');
    end

    C_in = uint8(C_in(:).');
    T_in = uint8(T_in(:).');

    if any(C_in ~= 0 & C_in ~= 1) || any(T_in ~= 0 & T_in ~= 1)
        error('C_in and T_in must be binary bit vectors containing only 0 and 1.');
    end

    BLOCKLEN = n/8;
    ell_len = ceil(log2(n));
    special_marker_pos = k + 1;

    num_blocks = length(C_in)/n;

    % ---------------------- outputs ----------------------
    c_decoded = -1 * ones(1, length(C_in), 'int8');
    M_out = uint8([]);
    valid = false;

    abandoned = zeros(1, num_blocks);
    n_guesses = zeros(1, num_blocks);

    info = struct();
    info.tag_n_guesses = 0;
    info.tag_abandoned = 0;
    info.tag_hat = [];
    info.tag_calc = [];
    info.decode_failure_reason = '';
    info.special_marker_pos = special_marker_pos;

    if taglen ~= n
        abandoned(:) = 1;
        info.decode_failure_reason = ...
            'This implementation supports tag GRAND decoding only when taglen = n.';
        return;
    end

    % ---------------------- initialize L and H ----------------------
    zero_bytes = zeros(1, BLOCKLEN, 'uint8');
    L = aes_enc_bytes(key_hex, zero_bytes);

    N_bytes = nonce_to_bytes(N, n);
    H = aes_enc_bytes(key_hex, bitxor(N_bytes, L));

    Z = rocb_offsets(L, H, num_blocks);

    Checksum = zeros(1, BLOCKLEN, 'uint8');

    info.L = L;
    info.H = H;
    info.Z = Z;
    info.ell_len = ell_len;
    info.num_blocks = num_blocks;

    % ================================================================
    % 1) Blocks 1,...,num_blocks-2: regular CRC codewords.
    % ================================================================
    if num_blocks > 1
        for i = 1:num_blocks-2
            C_block = get_block(C_in, i, n);

            [decoded, D_bits, C_hat, noise_found, guesses, block_kind] = ...
                rocb_grand_block_crc( ...
                    key_hex, C_block, L, Z(i,:), ...
                    'regular', k, ell_len, max_query, H_crc, []);

            n_guesses(i) = guesses;
            c_decoded((i-1)*n+1:i*n) = int8(C_hat);
            info.block_noise_found{i} = noise_found; %#ok<AGROW>
            info.block_kind{i} = block_kind; %#ok<AGROW>

            if decoded == 0
                abandoned(i) = 1;
                info.decode_failure_reason = ...
                    sprintf('GRAND abandoned regular CRC block %d.', i);
                M_out = [];
                return;
            end

            M_out = [M_out, D_bits(1:k)]; %#ok<AGROW>
            Checksum = bitxor(Checksum, bits_to_bytes(D_bits));
        end

        % ============================================================
        % 2) Penultimate block:
        %    either regular CRC or first split special CRC block.
        % ============================================================
        i = num_blocks - 1;
        C_block = get_block(C_in, i, n);

        [decoded, B, C_hat, noise_found, guesses, block_kind] = ...
            rocb_grand_block_crc( ...
                key_hex, C_block, L, Z(i,:), ...
                'regular_or_special', k, ell_len, ...
                max_query, H_crc, []);

        n_guesses(i) = guesses;
        c_decoded((i-1)*n+1:i*n) = int8(C_hat);
        info.block_noise_found{i} = noise_found;
        info.block_kind{i} = block_kind;

        if decoded == 0
            abandoned(i) = 1;
            info.decode_failure_reason = ...
                'GRAND abandoned penultimate CRC block.';
            M_out = [];
            return;
        end

        if strcmp(block_kind,'regular')
            M_out = [M_out, B(1:k)]; %#ok<AGROW>
            info.penultimate_type = 'regular_crc';

        elseif strcmp(block_kind,'special')
            U = B(1:k);

            len_bits = U(k-ell_len+1:k);
            n_rand = bits_to_uint(len_bits);

            % This has already been checked by the GRAND membership test.
            msg_len = k - ell_len - n_rand;

            M_out = [M_out, U(1:msg_len)]; %#ok<AGROW>
            info.penultimate_type = 'split_first_crc_special';
            info.penultimate_n_rand = n_rand;

        else
            error('Internal error: unexpected penultimate block kind.');
        end

        Checksum = bitxor(Checksum, bits_to_bytes(B));
    end

    % ================================================================
    % 3) Final ciphertext block: special CRC codeword only.
    % ================================================================
    i = num_blocks;
    C_block = get_block(C_in, i, n);

    [decoded, B_last, C_hat, noise_found, guesses, block_kind] = ...
        rocb_grand_block_crc( ...
            key_hex, C_block, L, Z(i,:), ...
            'special', k, ell_len, max_query, H_crc, []);

    n_guesses(i) = guesses;
    c_decoded((i-1)*n+1:i*n) = int8(C_hat);
    info.block_noise_found{i} = noise_found;
    info.block_kind{i} = block_kind;

    if decoded == 0
        abandoned(i) = 1;
        info.decode_failure_reason = 'GRAND abandoned final CRC block.';
        M_out = [];
        return;
    end

    U_last = B_last(1:k);
    len_bits = U_last(k-ell_len+1:k);
    n_rand_last = bits_to_uint(len_bits);

    % This has already been checked in the membership test.
    msg_len_last = k - ell_len - n_rand_last;
    M_last = U_last(1:msg_len_last);

    M_out = [M_out, M_last];
    Checksum = bitxor(Checksum, bits_to_bytes(B_last));

    info.final_n_rand = n_rand_last;
    info.final_B_last = B_last;
    info.Checksum = Checksum;

    % ================================================================
    % 4) Decode and verify the R-OCB tag.
    %
    % This part remains the same as zero-padding R-OCB.  The tag is not
    % replaced by a CRC codeword.
    % ================================================================
    Checksum_bits = bytes_to_bits(Checksum);
    P_tag = Checksum_bits(end-(n-k)+1:end);

    [tag_decoded, T_hat, tag_C_hat, tag_noise_found, tag_guesses] = ...
        rocb_grand_block_crc( ...
            key_hex, T_in, L, Z(num_blocks,:), ...
            'tag', k, ell_len, max_query, H_crc, P_tag);

    info.tag_n_guesses = tag_guesses;
    info.tag_noise_found = tag_noise_found;
    info.tag_hat = T_hat;
    info.tag_C_hat = tag_C_hat;
    info.P_tag = P_tag;

    if tag_decoded == 0
        info.tag_abandoned = 1;
        info.decode_failure_reason = 'GRAND abandoned tag.';
        M_out = [];
        return;
    end

    T_calc_bytes = aes_enc_bytes( ...
        key_hex, bitxor(Checksum, Z(num_blocks,:)));

    T_calc_bits = bytes_to_bits(T_calc_bytes);
    info.tag_calc = T_calc_bits;

    valid = isequal(uint8(T_hat), uint8(T_calc_bits(1:taglen)));

    if ~valid
        info.decode_failure_reason = ...
            'Decoded tag does not match calculated tag.';
        M_out = [];
    end
end

% =====================================================================
% GRAND helper
% =====================================================================
function [decoded, D_or_T_hat, C_hat, noise_found, n_guesses, block_kind] = ...
    rocb_grand_block_crc( ...
        key_hex, C_block, L, Z, membership, ...
        k, ell_len, max_query, H_crc, P_tag)

    n = length(C_block);

    decoded = 0;
    D_or_T_hat = [];
    C_hat = -1 * ones(1,n,'int8');
    noise_found = NaN(1,n);
    n_guesses = 0;
    block_kind = '';

    [err_loc_vec, err_vec, ~] = gen_next_err(n);
    L_inv = gf128_mul_x_inv(L);

    while decoded == 0 && n_guesses < max_query
        n_guesses = n_guesses + 1;

        C_hat = int8(mod(double(C_block) + double(err_vec), 2));
        C_hat_bytes = bits_to_bytes(uint8(C_hat));

        try
            if strcmp(membership,'tag')
                D_bytes = aes_dec_bytes(key_hex, C_hat_bytes);
                D_bytes = bitxor(D_bytes, Z);

            elseif strcmp(membership,'special')
                tmp = bitxor(C_hat_bytes, Z);
                D_bytes = aes_dec_bytes(key_hex, tmp);
                D_bytes = bitxor(D_bytes, L_inv);
                D_bytes = bitxor(D_bytes, Z);

            else
                % regular or regular_or_special use ordinary R-OCB inverse.
                tmp = bitxor(C_hat_bytes, Z);
                D_bytes = aes_dec_bytes(key_hex, tmp);
                D_bytes = bitxor(D_bytes, Z);
            end
        catch
            [err_loc_vec, err_vec, weight] = ...
                gen_next_err(n, err_loc_vec);

            if isinf(weight)
                break;
            end
            continue;
        end

        D_bits = bytes_to_bits(D_bytes);

        % -------------------------------------------------------------
        % Membership test
        % -------------------------------------------------------------
        if strcmp(membership,'regular')
            if crc_regular_is_valid(D_bits, H_crc)
                decoded = 1;
                block_kind = 'regular';
            end

        elseif strcmp(membership,'special')
            if crc_special_is_valid(D_bits, H_crc, k, ell_len)
                decoded = 1;
                block_kind = 'special';
            end

        elseif strcmp(membership,'regular_or_special')
            if crc_regular_is_valid(D_bits, H_crc)
                decoded = 1;
                block_kind = 'regular';

            elseif crc_special_is_valid(D_bits, H_crc, k, ell_len)
                decoded = 1;
                block_kind = 'special';
            end

        elseif strcmp(membership,'tag')
            pad_len = length(P_tag);
            P_test = D_bits(end-pad_len+1:end);

            if isequal(uint8(P_test), uint8(P_tag))
                decoded = 1;
                block_kind = 'tag';
            end

        else
            error('Unknown CRC-R-OCB membership type: %s.', membership);
        end

        if decoded
            noise_found = uint8(err_vec);

            if strcmp(membership,'tag')
                D_or_T_hat = uint8(C_hat);
            else
                D_or_T_hat = D_bits;
            end
            return;
        end

        [err_loc_vec, err_vec, weight] = gen_next_err(n, err_loc_vec);

        if isinf(weight)
            break;
        end
    end
end

% =====================================================================
% CRC membership tests
% =====================================================================
function ok = crc_regular_is_valid(D_bits, H_crc)
    syndrome = mod(double(H_crc) * double(D_bits(:)), 2);
    ok = all(syndrome == 0);
end

function ok = crc_special_is_valid(D_bits, H_crc, k, ell_len)
    % Undo the special coset marker before the CRC syndrome test.
    B_crc = uint8(D_bits(:).');
    marker_pos = k + 1;
    B_crc(marker_pos) = uint8(1 - B_crc(marker_pos));

    syndrome = mod(double(H_crc) * double(B_crc(:)), 2);

    if ~all(syndrome == 0)
        ok = false;
        return;
    end

    % Also enforce a structurally valid special-block length field.
    U = uint8(D_bits(1:k));
    len_bits = U(k-ell_len+1:k);
    n_rand = bits_to_uint(len_bits);

    ok = (n_rand <= k-ell_len);
end

% =====================================================================
% Local helpers
% =====================================================================
function block = get_block(bits, i, n)
    block = bits((i-1)*n+1:i*n);
end

function out_bytes = aes_enc_bytes(key_hex, in_bytes)
    in_hex = upper(reshape(dec2hex(uint8(in_bytes), 2).', 1, []));
    out_hex = Cipher(key_hex, in_hex);
    out_bytes = uint8(hex2dec(reshape(out_hex, 2, []).')).';
end

function out_bytes = aes_dec_bytes(key_hex, in_bytes)
    in_hex = upper(reshape(dec2hex(uint8(in_bytes), 2).', 1, []));
    out_hex = InvCipher(key_hex, in_hex);
    out_bytes = uint8(hex2dec(reshape(out_hex, 2, []).')).';
end

function bytes = bits_to_bytes(bits)
    bits = uint8(bits(:).');

    if mod(length(bits),8) ~= 0
        error('Bit vector length must be divisible by 8.');
    end

    bytes = uint8(bi2de(reshape(bits, 8, []).', 'left-msb')).';
end

function bits = bytes_to_bits(bytes)
    bytes = uint8(bytes(:).');
    bits = uint8(reshape(de2bi(bytes, 8, 'left-msb').', 1, []));
end

function value = bits_to_uint(bits)
    bits = uint8(bits(:).');
    value = 0;

    for idx = 1:length(bits)
        value = 2*value + double(bits(idx));
    end
end

function N_bytes = nonce_to_bytes(N, n)
    if ischar(N) || isstring(N)
        N_hex = char(N);

        if length(N_hex) ~= n/4
            error('Hex nonce must contain n/4 hex characters.');
        end

        N_bytes = uint8(hex2dec(reshape(N_hex, 2, []).')).';
    else
        N = uint8(N(:).');

        if length(N) ~= n || any(N ~= 0 & N ~= 1)
            error('Bit-vector nonce must contain exactly n binary bits.');
        end

        N_bytes = bits_to_bytes(N);
    end
end

function Z = rocb_offsets(L, H, num_offsets)
    if exist('compute_rocb_offsets', 'file') == 2
        Z_existing = compute_rocb_offsets(L, H, num_offsets);
        Z = uint8(Z_existing);
        return;
    end

    Z = zeros(num_offsets, length(L), 'uint8');

    for i = 1:num_offsets
        gammaL = rocb_gamma_times_L(L, i);
        Z(i,:) = bitxor(gammaL, H);
    end
end

function gammaL = rocb_gamma_times_L(L, i)
    gamma = bitxor(uint64(i), bitshift(uint64(i), -1));
    gammaL = gf128_mul_uint64(L, gamma);
end

function out = gf128_mul_uint64(block, scalar)
    out = zeros(size(block), 'uint8');
    term = uint8(block);

    while scalar > 0
        if bitand(scalar, uint64(1))
            out = bitxor(out, term);
        end

        scalar = bitshift(scalar, -1);

        if scalar > 0
            term = gf128_mul_x(term);
        end
    end
end

function y = gf128_mul_x(x)
    x = uint8(x(:).');
    y = zeros(size(x), 'uint8');
    carry = uint8(0);

    for j = length(x):-1:1
        new_carry = bitshift(bitand(x(j), uint8(128)), -7);
        y(j) = bitor(bitshift(x(j), 1), carry);
        carry = new_carry;
    end

    if carry
        y(end) = bitxor(y(end), uint8(hex2dec('87')));
    end
end

function x = gf128_mul_x_inv(y)
    y = uint8(y(:).');

    if bitand(y(end), uint8(1)) == 0
        x = right_shift_one(y);
    else
        y2 = y;
        y2(end) = bitxor(y2(end), uint8(hex2dec('87')));
        x = right_shift_one(y2);
        x(1) = bitor(x(1), uint8(128));
    end
end

function y = right_shift_one(x)
    x = uint8(x(:).');
    y = zeros(size(x), 'uint8');
    carry = uint8(0);

    for j = 1:length(x)
        new_carry = bitand(x(j), uint8(1));
        y(j) = bitor(bitshift(x(j), -1), bitshift(carry, 7));
        carry = new_carry;
    end
end

% =====================================================================
% GRAND hard-decision error-pattern generator
% =====================================================================
function [err_loc_vec, err_vec, weight] = gen_next_err(n, err_loc_vec)
    mask_vec = 1:n;
    err_vec = zeros(1,n,'uint8');

    if ~exist('err_loc_vec','var')
        err_loc_vec = [];
        weight = 0;
        return;
    end

    err_loc_vec = increase_error(mask_vec, err_loc_vec);

    if isempty(err_loc_vec)
        weight = inf;
        return;
    end

    weight = length(err_loc_vec);
    err_vec(err_loc_vec) = 1;
end

function err_loc_vec = increase_error(mask_vec, err_loc_vec)
    max_possible_err_loc = length(mask_vec);
    success = false;

    for ii = length(err_loc_vec):-1:1
        new_err_loc_vec = err_loc_vec;

        if new_err_loc_vec(ii) == max(mask_vec)
            continue;
        end

        success = true;
        ind_in_mask = find(mask_vec == new_err_loc_vec(ii), 1);
        new_err_loc_vec(ii) = mask_vec(ind_in_mask+1);

        for jj = ii+1:length(err_loc_vec)
            new_ind = ind_in_mask + jj - ii + 1;

            if new_ind <= max_possible_err_loc
                new_err_loc_vec(jj) = mask_vec(new_ind);
            else
                success = false;
                break;
            end
        end

        if success
            err_loc_vec = new_err_loc_vec;
            return;
        end
    end

    if ~success
        if length(err_loc_vec) == max_possible_err_loc
            err_loc_vec = [];
            return;
        end

        err_loc_vec = mask_vec(1:length(err_loc_vec)+1);
    end
end
