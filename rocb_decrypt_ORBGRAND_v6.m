function [c_decoded, M_out, valid, abandoned, n_guesses, info] = rocb_decrypt_ORBGRAND_v6(key_hex, N, C_in, T_in, taglen, n, k, max_query, P)
%ROCB_DECRYPT_ORBGRAND_V6  Joint ROCB decryption-decoding using ORBGRAND.
%
%   [c_decoded,M_out,valid,abandoned,n_guesses,info] = ...
%       rocb_decrypt_ORBGRAND_v6(key_hex,N,C_in,T_in,taglen,n,k,max_query,P)
%
% This implementation follows the ROCB decryption pseudocode.  It is the
% ORBGRAND/AWGN counterpart of rocb_decrypt_GRAND_v6: the ROCB inverse
% operations, padding tests, checksum update, and tag verification are the
% same, but the candidate error patterns are ordered according to bit
% reliabilities derived from the soft channel output.
%
% Inputs
%   key_hex   : AES key as hex string, used by Cipher/InvCipher.
%   N         : nonce, either n-bit vector or n/4-character hex string.
%   C_in      : noisy ciphertext, without the tag. For AWGN/BPSK this is a
%               soft vector. Hard decision uses bit = (soft < 0). A binary
%               vector is also accepted.
%   T_in      : noisy tag. The pseudocode supports tag decoding only when
%               taglen = n. T_in may be soft or binary.
%   taglen    : tag length in bits. Must be n for this pseudocode version.
%   n         : block length in bits, e.g., 128.
%   k         : information length in bits, e.g., 112, 116, 120.
%   max_query : maximum ORBGRAND guesses per block.
%   P         : fixed redundancy padding of length n-k bits. Optional;
%               default is all zeros.
%
% Outputs
%   c_decoded : decoded/corrected ciphertext bits, without the tag.
%   M_out     : decoded message bits, or [] if invalid/abandoned.
%   valid     : true if all blocks decoded and the tag verified.
%   abandoned : per-ciphertext-block abandonment indicators.
%   n_guesses : per-ciphertext-block ORBGRAND query counts.
%   info      : debug information, including tag ORBGRAND statistics.

    % ---------------------- input checks ----------------------
    if nargin < 9 || isempty(P)
        P = zeros(1, n-k, 'uint8');
    end

    validateattributes(n, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(k, {'numeric'}, {'scalar','integer','positive','<',n});
    validateattributes(taglen, {'numeric'}, {'scalar','integer','positive','<=',n});
    % validateattributes(max_query, {'numeric'}, {'scalar','integer','positive'});

    if mod(n,8) ~= 0
        error('n must be a multiple of 8 because AES/Cipher is byte based.');
    end
    if length(P) ~= n-k
        error('P must contain exactly n-k bits.');
    end
    if isempty(C_in) || mod(length(C_in), n) ~= 0
        error('C_in length must be a positive multiple of n.');
    end
    if length(T_in) ~= taglen
        error('T_in must contain exactly taglen entries.');
    end

    C_in = double(C_in(:).');
    T_in = double(T_in(:).');
    P = uint8(P(:).');
    if any(P ~= 0 & P ~= 1)
        error('P must be a binary bit vector containing only 0 and 1.');
    end

    BLOCKLEN = n/8;
    ell_len = ceil(log2(n));        % 7 when n = 128
    if ell_len > k
        error('k must be at least ceil(log2(n)).');
    end
    if n-k < 1
        error('The construction requires at least one padding bit, i.e., k < n.');
    end

    num_blocks = length(C_in) / n;  % \bar m_C in the pseudocode

    % Default outputs. They are overwritten as decoding succeeds.
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

    % The pseudocode only defines tag GRAND/ORBGRAND decoding for tau = n.
    if taglen ~= n
        abandoned(:) = 1;
        info.decode_failure_reason = 'The pseudocode supports tag decoding only when taglen = n.';
        return;
    end

    % ---------------------- initialize L and H ----------------------
    zero_bytes = zeros(1, BLOCKLEN, 'uint8');
    L = aes_enc_bytes(key_hex, zero_bytes);          % L <- E_K(0^n)

    N_bytes = nonce_to_bytes(N, n);
    H = aes_enc_bytes(key_hex, bitxor(N_bytes, L));  % H <- E_K(N xor L)

    % ---------------------- offsets and padding ----------------------
    Z = rocb_offsets(L, H, num_blocks);

    % P_tilde <- complement of first padding bit, followed by the rest.
    % For the usual P = 0^{n-k}, this is 1 || 0^{n-k-1}.
    P_tilde = [uint8(1 - P(1)), P(2:end)];

    Checksum = zeros(1, BLOCKLEN, 'uint8');

    info.L = L;
    info.H = H;
    info.Z = Z;
    info.P_tilde = P_tilde;
    info.ell_len = ell_len;
    info.num_blocks = num_blocks;

    % ================================================================
    % 1) Decode blocks 1,...,num_blocks-2 as full/regular blocks with P.
    % ================================================================
    if num_blocks > 1
        for i = 1:num_blocks-2
            C_block_soft = get_block(C_in, i, n);

            [decoded, D_bits, C_hat, noise_found, guesses] = ...
                rocb_orbgrand_block(key_hex, C_block_soft, L, Z(i,:), {P}, 'full', k, max_query);

            n_guesses(i) = guesses;
            c_decoded((i-1)*n+1:i*n) = int8(C_hat);
            info.block_noise_found{i} = noise_found; %#ok<AGROW>

            if decoded == 0
                abandoned(i) = 1;
                info.decode_failure_reason = sprintf('ORBGRAND abandoned full block %d.', i);
                valid = false;
                M_out = [];
                return;
            end

            M_out = [M_out, D_bits(1:k)]; %#ok<AGROW>
            Checksum = bitxor(Checksum, bits_to_bytes(D_bits));
        end

        % ============================================================
        % 2) Decode the penultimate block. It may be either:
        %    - a full block ending with P, or
        %    - the first part of a split last block ending with P_tilde.
        % ============================================================
        i = num_blocks - 1;
        C_block_soft = get_block(C_in, i, n);

        [decoded, B, C_hat, noise_found, guesses] = ...
            rocb_orbgrand_block(key_hex, C_block_soft, L, Z(i,:), {P, P_tilde}, 'full', k, max_query);

        n_guesses(i) = guesses;
        c_decoded((i-1)*n+1:i*n) = int8(C_hat);
        info.block_noise_found{i} = noise_found;

        if decoded == 0
            abandoned(i) = 1;
            info.decode_failure_reason = 'ORBGRAND abandoned penultimate block.';
            valid = false;
            M_out = [];
            return;
        end

        P_pen = B(end-(n-k)+1:end);
        if isequal(P_pen, P)
            % B = M || P
            M_out = [M_out, B(1:k)]; %#ok<AGROW>
            info.penultimate_type = 'full_P';
        elseif isequal(P_pen, P_tilde)
            % B = M_part || P_rand || n_P_rand || P_tilde.
            len_bits = B(k-ell_len+1:k);
            n_rand = bits_to_uint(len_bits);
            if n_rand > k-ell_len
                abandoned(i) = 1;
                info.decode_failure_reason = 'Invalid random-padding length in penultimate block.';
                valid = false;
                M_out = [];
                return;
            end
            msg_len = k - ell_len - n_rand;
            M_out = [M_out, B(1:msg_len)]; %#ok<AGROW>
            info.penultimate_type = 'split_first_P_tilde';
            info.penultimate_n_rand = n_rand;
        else
            abandoned(i) = 1;
            info.decode_failure_reason = 'Penultimate block padding is neither P nor P_tilde.';
            valid = false;
            M_out = [];
            return;
        end

        Checksum = bitxor(Checksum, bits_to_bytes(B));
    end

    % ================================================================
    % 3) Decode the final ciphertext block with the special last rule.
    %    B_last = M_last || P_rand || n_P_rand || P_tilde.
    % ================================================================
    i = num_blocks;
    C_block_soft = get_block(C_in, i, n);

    [decoded, B_last, C_hat, noise_found, guesses] = ...
        rocb_orbgrand_block(key_hex, C_block_soft, L, Z(i,:), {P_tilde}, 'last', k, max_query);

    n_guesses(i) = guesses;
    c_decoded((i-1)*n+1:i*n) = int8(C_hat);
    info.block_noise_found{i} = noise_found;

    if decoded == 0
        abandoned(i) = 1;
        info.decode_failure_reason = 'ORBGRAND abandoned final block.';
        valid = false;
        M_out = [];
        return;
    end

    len_bits = B_last(k-ell_len+1:k);
    n_rand_last = bits_to_uint(len_bits);
    if n_rand_last > k-ell_len
        abandoned(i) = 1;
        info.decode_failure_reason = 'Invalid random-padding length in final block.';
        valid = false;
        M_out = [];
        return;
    end

    msg_len_last = k - ell_len - n_rand_last;
    M_last = B_last(1:msg_len_last);
    M_out = [M_out, M_last];
    Checksum = bitxor(Checksum, bits_to_bytes(B_last));

    info.final_n_rand = n_rand_last;
    info.final_B_last = B_last;
    info.Checksum = Checksum;

    % ================================================================
    % 4) Decode and verify the tag.
    % ================================================================
    Checksum_bits = bytes_to_bits(Checksum);
    P_tag = Checksum_bits(end-(n-k)+1:end);

    [tag_decoded, T_hat, tag_C_hat, tag_noise_found, tag_guesses] = ...
        rocb_orbgrand_block(key_hex, T_in, L, Z(num_blocks,:), {P_tag}, 'tag', k, max_query);

    info.tag_n_guesses = tag_guesses;
    info.tag_noise_found = tag_noise_found;
    info.tag_hat = T_hat;
    info.tag_C_hat = tag_C_hat;
    info.P_tag = P_tag;

    if tag_decoded == 0
        info.tag_abandoned = 1;
        info.decode_failure_reason = 'ORBGRAND abandoned tag.';
        valid = false;
        M_out = [];
        % return;
    end

    T_calc_bytes = aes_enc_bytes(key_hex, bitxor(Checksum, Z(num_blocks,:)));
    T_calc_bits = bytes_to_bits(T_calc_bytes);
    info.tag_calc = T_calc_bits;

    valid = isequal(uint8(T_hat), uint8(T_calc_bits(1:taglen)));
    if ~valid
        info.decode_failure_reason = 'Decoded tag does not match calculated tag.';
        M_out = [];
    end
end

% =====================================================================
% ROCB.ORBGRAND helper
% =====================================================================
function [decoded, D_or_T_hat, C_hat, noise_found, n_guesses] = ...
    rocb_orbgrand_block(key_hex, C_block_soft, L, Z, P_acc, type, k, max_query)
%ROCB_ORBGRAND_BLOCK  ORBGRAND membership test for one ROCB block.
%
% type = 'full': D = InvCipher(C_hat xor Z) xor Z
% type = 'last': D = InvCipher(C_hat xor Z) xor L*x^{-1} xor Z
% type = 'tag' : D = InvCipher(C_hat) xor Z; return C_hat on success
%
% ORBGRAND difference from GRAND:
%   - C_block_soft is hard-demodulated to a binary reference word.
%   - Candidate error locations are generated in reliability order.
%   - The ROCB inverse operation and padding membership check are unchanged.

    n = length(C_block_soft);
    decoded = 0;
    D_or_T_hat = [];
    C_hat = -1 * ones(1, n, 'int8');
    noise_found = NaN(1, n);
    n_guesses = 0;

    C_ref = hard_demodulate(C_block_soft);
    reliability = bit_reliability(C_block_soft);
    [~, reliability_order] = sort(reliability, 'ascend');  % least reliable first

    L_inv = gf128_mul_x_inv(L);

    % ---------------- first query: zero-noise candidate ----------------
    [decoded, D_or_T_hat, C_hat, noise_found] = try_candidate( ...
        key_hex, C_ref, zeros(1,n,'uint8'), L, L_inv, Z, P_acc, type, k);
    n_guesses = 1;
    if decoded == 1 || n_guesses >= max_query
        return;
    end

    % ---------------- ORBGRAND reliability-ordered queries -------------
    W = 0;
    while decoded == 0 && n_guesses < max_query
        W = W + 1;
        w = ceil((1 + 2*n - sqrt((1 + 2*n)^2 - 8*W)) / 2);

        while decoded == 0 && n_guesses < max_query && ...
              w <= floor((sqrt(1 + 8*W) - 1) / 2)

            noise_locations = landslide(W, w, n);

            for jj = 1:size(noise_locations, 1)
                if n_guesses >= max_query
                    break;
                end

                err_vec = zeros(1, n, 'uint8');
                ordered_locations = noise_locations(jj, :);
                ordered_locations = ordered_locations(ordered_locations >= 1 & ordered_locations <= n);

                % Map ORBGRAND ordered-domain locations to original bit positions.
                err_vec(reliability_order(ordered_locations)) = 1;

                [decoded, D_or_T_hat, C_hat, noise_found] = try_candidate( ...
                    key_hex, C_ref, err_vec, L, L_inv, Z, P_acc, type, k);
                n_guesses = n_guesses + 1;

                if decoded == 1
                    return;
                end
            end

            w = w + 1;
        end
    end
end

function [decoded, D_or_T_hat, C_hat, noise_found] = try_candidate( ...
    key_hex, C_ref, err_vec, L, L_inv, Z, P_acc, type, k)

    n = length(C_ref);
    decoded = 0;
    D_or_T_hat = [];
    C_hat = int8(mod(double(C_ref) + double(err_vec), 2));
    noise_found = NaN(1, n);

    try
        C_hat_bytes = bits_to_bytes(uint8(C_hat));

        switch type
            case 'full'
                tmp = bitxor(C_hat_bytes, Z);
                D_bytes = aes_dec_bytes(key_hex, tmp);
                D_bytes = bitxor(D_bytes, Z);

            case 'last'
                tmp = bitxor(C_hat_bytes, Z);
                D_bytes = aes_dec_bytes(key_hex, tmp);
                D_bytes = bitxor(D_bytes, L_inv);
                D_bytes = bitxor(D_bytes, Z);

            case 'tag'
                D_bytes = aes_dec_bytes(key_hex, C_hat_bytes);
                D_bytes = bitxor(D_bytes, Z);

            otherwise
                error('Unknown ROCB ORBGRAND block type: %s.', type);
        end
    catch
        return;
    end

    D_bits = bytes_to_bits(D_bytes);
    pad_len = length(P_acc{1});
    P_test = D_bits(end-pad_len+1:end);

    if padding_is_acceptable(P_test, P_acc)
        decoded = 1;
        noise_found = uint8(err_vec);
        if strcmp(type, 'tag')
            % The pseudocode returns the corrected tag block C_hat.
            D_or_T_hat = uint8(C_hat);
        else
            D_or_T_hat = D_bits;
        end
    end
end

function bits = hard_demodulate(x)
    x = double(x(:).');
    if all(x == 0 | x == 1)
        bits = uint8(x);
    else
        bits = uint8(x < 0);
    end
end

function r = bit_reliability(x)
    x = double(x(:).');
    if all(x == 0 | x == 1)
        % No soft information is available. Keep the natural bit order.
        r = ones(size(x));
    else
        r = abs(x);
    end
end

function ok = padding_is_acceptable(P_test, P_acc)
    ok = false;
    for j = 1:numel(P_acc)
        if isequal(uint8(P_test), uint8(P_acc{j}))
            ok = true;
            return;
        end
    end
end

% =====================================================================
% Local helper functions
% =====================================================================
function block = get_block(x, i, n)
    block = x((i-1)*n+1:i*n);
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
    if mod(length(bits), 8) ~= 0
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
    % Prefer the existing project helper, if it exists, so encryption and
    % decryption remain compatible with each other.
    if exist('compute_rocb_offsets', 'file') == 2
        Z_existing = compute_rocb_offsets(L, H, num_offsets);
        Z = uint8(Z_existing);
        return;
    end

    % Fallback: gamma_i * L xor H, using Gray-code gamma_i.
    Z = zeros(num_offsets, length(L), 'uint8');
    for i = 1:num_offsets
        gammaL = rocb_gamma_times_L(L, i);
        Z(i,:) = bitxor(gammaL, H);
    end
end

function gammaL = rocb_gamma_times_L(L, i)
    gamma = bitxor(uint64(i), bitshift(uint64(i), -1)); % Gray-code gamma_i
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
