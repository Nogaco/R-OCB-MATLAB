function [C_out, T_out, info] = rocb_encrypt_CRC_v1(key_hex, N, M, taglen, n, k, CRC)
%ROCB_ENCRYPT_CRC_V1  R-OCB encryption with CRC redundancy.
%
%   [C_out,T_out,info] = ...
%       rocb_encrypt_CRC_v1(key_hex,N,M,taglen,n,k,CRC)

% Inputs
%   key_hex : AES key as hex string, used by Cipher(key_hex, block_hex)
%   N       : nonce, either n-bit vector or n/4-character hex string
%   M       : row vector of message bits, values 0/1
%   taglen  : tag length in bits
%   n       : block length in bits, e.g., 128
%   k       : information length in bits, e.g., 112, 116, 120
%
% Outputs
%   C_out : ciphertext bits, without the tag
%   T_out : first taglen bits of the authentication tag
%   info  : struct with useful debug information
%
% CRC must contain:
%   CRC.G : k-by-n systematic CRC generator matrix
%   CRC.H : parity-check matrix (not used by the encryptor, but kept in CRC)
%
% Regular R-OCB block:
%   U = k payload bits
%   B = CRC_encode(U) = [U || CRC(U)]
%
% Special/split-last R-OCB block:
%   U = [message-part || random-padding || length-field]  (k bits)
%   B = CRC_encode(U)
%   Then bit k+1 is complemented to mark the special block.
%
% The AES/R-OCB transformations, offsets, checksum, and tag construction
% are unchanged from the zero-padding R-OCB implementation.

    % ---------------------- input checks ----------------------
    validateattributes(n, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(k, {'numeric'}, {'scalar','integer','positive','<',n});
    validateattributes(taglen, {'numeric'}, {'scalar','integer','positive','<=',n});

    if mod(n,8) ~= 0
        error('n must be a multiple of 8 because AES/Cipher is byte based.');
    end

    if ~isstruct(CRC) || ~isfield(CRC,'G')
        error('CRC must be a struct containing CRC.G.');
    end

    G_crc = double(CRC.G);

    if size(G_crc,1) ~= k || size(G_crc,2) ~= n
        error('CRC.G must have size k-by-n = %d-by-%d.', k, n);
    end

    M = uint8(M(:).');
    if any(M ~= 0 & M ~= 1)
        error('M must be a binary bit vector containing only 0 and 1.');
    end

    BLOCKLEN = n/8;
    ell_len = ceil(log2(n));        % 7 when n = 128

    if ell_len > k
        error('k must be at least ceil(log2(n)).');
    end

    % ---------------------- initialize L and H ----------------------
    zero_bytes = zeros(1, BLOCKLEN, 'uint8');
    L = aes_enc_bytes(key_hex, zero_bytes);          % L <- E_K(0^n)

    N_bytes = nonce_to_bytes(N, n);
    H = aes_enc_bytes(key_hex, bitxor(N_bytes, L));  % H <- E_K(N xor L)

    % ---------------------- partition M ----------------------
    mlen_bits = length(M);
    mbar = max(1, ceil(mlen_bits/k));

    % Regular blocks are M[1],...,M[mbar-1].
    regular_count = mbar - 1;

    Checksum = zeros(1, BLOCKLEN, 'uint8');

    % Determine the final message part.
    last_start = regular_count*k + 1;
    if last_start <= mlen_bits
        M_last = M(last_start:end);
    else
        M_last = uint8([]);
    end

    n_last = length(M_last);
    split_last = (n_last > k - ell_len);

    % The final pseudocode branch creates either one or two ciphertext blocks.
    num_cipher_blocks = mbar + split_last;

    Z = rocb_offsets(L, H, num_cipher_blocks);

    C_out = zeros(1, num_cipher_blocks*n, 'uint8');
    C_bits_idx = 1;

    % ================================================================
    % 1) Regular blocks: U -> systematic CRC codeword B
    % ================================================================
    for i = 1:regular_count
        U = M((i-1)*k + 1 : i*k);

        B = crc_encode_payload(U, G_crc, k, n);

        [C_block_bits, B_bytes] = ...
            encrypt_rocb_block(key_hex, B, Z(i,:));

        C_out(C_bits_idx:C_bits_idx+n-1) = C_block_bits;
        C_bits_idx = C_bits_idx + n;

        Checksum = bitxor(Checksum, B_bytes);
    end

    % ================================================================
    % 2) Special final block(s)
    %
    % Instead of fixed P_tilde, CRC-R-OCB uses the coset marker obtained
    % by complementing the first redundancy bit, position k+1.
    % ================================================================
    L_x_inv = gf128_mul_x_inv(L);
    special_marker_pos = k + 1;

    if ~split_last
        % ------------------------------------------------------------
        % One special final block.
        % U = M_last || P_rand || len(n_rand), exactly k bits.
        % ------------------------------------------------------------
        n_rand = k - ell_len - n_last;
        P_rand = uint8(randi([0,1], 1, n_rand));
        len_bits = uint_to_bits(n_rand, ell_len);

        U = [M_last, P_rand, len_bits];

        if length(U) ~= k
            error('Internal error: final CRC payload is not k bits.');
        end

        B = crc_encode_payload(U, G_crc, k, n);
        B(special_marker_pos) = uint8(1 - B(special_marker_pos));

        Z_last_masked = bitxor(Z(mbar,:), L_x_inv);

        [C_block_bits, B_bytes] = ...
            encrypt_rocb_last_block( ...
                key_hex, B, Z_last_masked, Z(mbar,:));

        C_out(C_bits_idx:C_bits_idx+n-1) = C_block_bits;

        Checksum = bitxor(Checksum, B_bytes);

        last_tag_offset_index = mbar;
        B_last = B;

    else
        % ------------------------------------------------------------
        % Split final message block into two encoded blocks.
        % First split block uses the ordinary R-OCB block transform but
        % is marked as special.
        % ------------------------------------------------------------
        n_next = n_last - (k - ell_len);

        M_first  = M_last(1:k-ell_len);
        M_second = M_last(k-ell_len+1:end);

        % First split block: length field = 0.
        zero_len_bits = zeros(1, ell_len, 'uint8');
        U_first = [M_first, zero_len_bits];

        if length(U_first) ~= k
            error('Internal error: first split CRC payload is not k bits.');
        end

        B_first = crc_encode_payload(U_first, G_crc, k, n);
        B_first(special_marker_pos) = ...
            uint8(1 - B_first(special_marker_pos));

        [C_block_bits, B_first_bytes] = ...
            encrypt_rocb_block(key_hex, B_first, Z(mbar,:));

        C_out(C_bits_idx:C_bits_idx+n-1) = C_block_bits;
        C_bits_idx = C_bits_idx + n;

        Checksum = bitxor(Checksum, B_first_bytes);

        % Second split block: remaining bits + random padding + length.
        n_rand2 = k - ell_len - n_next;
        P_rand2 = uint8(randi([0,1], 1, n_rand2));
        len_bits2 = uint_to_bits(n_rand2, ell_len);

        U_last = [M_second, P_rand2, len_bits2];

        if length(U_last) ~= k
            error('Internal error: second split CRC payload is not k bits.');
        end

        B_last = crc_encode_payload(U_last, G_crc, k, n);
        B_last(special_marker_pos) = ...
            uint8(1 - B_last(special_marker_pos));

        Z_last_masked = bitxor(Z(mbar+1,:), L_x_inv);

        [C_block_bits, B_last_bytes] = ...
            encrypt_rocb_last_block( ...
                key_hex, B_last, Z_last_masked, Z(mbar+1,:));

        C_out(C_bits_idx:C_bits_idx+n-1) = C_block_bits;

        Checksum = bitxor(Checksum, B_last_bytes);

        last_tag_offset_index = mbar + 1;
    end

    % ---------------------- tag ----------------------
    tag_input = bitxor(Checksum, Z(last_tag_offset_index,:));
    T_bytes = aes_enc_bytes(key_hex, tag_input);
    T_bits = bytes_to_bits(T_bytes);
    T_out = T_bits(1:taglen);

    % ---------------------- debug info ----------------------
    if nargout > 2
        info = struct();
        info.L = L;
        info.H = H;
        info.Z = Z;
        info.mbar = mbar;
        info.regular_count = regular_count;
        info.split_last = split_last;
        info.n_last = n_last;
        info.ell_len = ell_len;
        info.special_marker_pos = special_marker_pos;
        info.B_last = B_last;
        info.Checksum = Checksum;
        info.last_tag_offset_index = last_tag_offset_index;
    end
end

% =====================================================================
% CRC helper
% =====================================================================
function B = crc_encode_payload(U, G_crc, k, n)
    U = uint8(U(:).');

    if length(U) ~= k
        error('CRC payload must contain exactly k bits.');
    end

    B = uint8(mod(double(U) * double(G_crc), 2));

    if length(B) ~= n
        error('CRC codeword must contain exactly n bits.');
    end

    % R-OCB relies on a systematic CRC generator:
    % the first k bits of the codeword are the payload U.
    if ~isequal(B(1:k), U)
        error(['CRC generator matrix is not systematic in the expected ' ...
               'form [I_k | parity].']);
    end
end

% =====================================================================
% R-OCB block encryption helpers
% =====================================================================
function [C_bits, B_bytes] = encrypt_rocb_block(key_hex, B_bits, Z_bytes)
    B_bits = uint8(B_bits(:).');
    Z_bytes = uint8(Z_bytes(:).');

    if mod(length(B_bits),8) ~= 0
        error('Internal block length must be a multiple of 8 bits.');
    end

    B_bytes = bits_to_bytes(B_bits);

    tmp_in = bitxor(B_bytes, Z_bytes);
    tmp_out = aes_enc_bytes(key_hex, tmp_in);
    C_bytes = bitxor(tmp_out, Z_bytes);

    C_bits = bytes_to_bits(C_bytes);
end

function [C_bits, B_bytes] = ...
    encrypt_rocb_last_block(key_hex, B_bits, Z_with_x_bytes, Z_bytes)

    B_bits = uint8(B_bits(:).');
    Z_bytes = uint8(Z_bytes(:).');
    Z_with_x_bytes = uint8(Z_with_x_bytes(:).');

    if mod(length(B_bits),8) ~= 0
        error('Internal block length must be a multiple of 8 bits.');
    end

    B_bytes = bits_to_bytes(B_bits);

    tmp_in = bitxor(B_bytes, Z_with_x_bytes);
    tmp_out = aes_enc_bytes(key_hex, tmp_in);
    C_bytes = bitxor(tmp_out, Z_bytes);

    C_bits = bytes_to_bits(C_bytes);
end

% =====================================================================
% Generic helpers
% =====================================================================
function out_bytes = aes_enc_bytes(key_hex, in_bytes)
    in_hex = upper(reshape(dec2hex(uint8(in_bytes), 2).', 1, []));
    out_hex = Cipher(key_hex, in_hex);
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

function bits = uint_to_bits(value, width)
    if value < 0 || value >= 2^width
        error('The value %d cannot be represented using %d bits.', value, width);
    end
    bits = uint8(dec2bin(value, width) - '0');
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
