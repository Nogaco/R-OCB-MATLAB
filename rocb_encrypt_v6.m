function [C_out, T_out, info] = rocb_encrypt_v6(key_hex, N, M, taglen, n, k, P)
%ROCB_ENCRYPT_V6  ROCB encryption according to the last-block pseudocode.
%
%   [C_out,T_out,info] = rocb_encrypt_v6(key_hex,N,M,taglen,n,k,P)
%
% Inputs
%   key_hex : AES key as hex string, used by Cipher(key_hex, block_hex)
%   N       : nonce, either n-bit vector or n/4-character hex string
%   M       : row vector of message bits, values 0/1
%   taglen  : tag length in bits
%   n       : block length in bits, e.g., 128
%   k       : information length in bits, e.g., 112, 116, 120
%   P       : fixed redundancy padding of length n-k bits. Optional;
%             default is all zeros.
%
% Outputs
%   C_out : ciphertext bits, without the tag
%   T_out : first taglen bits of the authentication tag
%   info  : struct with useful debug information
%
% Notes
%   1. This implementation treats the last message block exactly as in the
%      pseudocode. Therefore, even when length(M) is a multiple of k, the
%      last k-bit block is split so that a length field can be encoded.
%   2. The offset computation first tries to call an existing
%      compute_rocb_offsets(L,H,num_offsets) helper, if it exists on the
%      MATLAB path. If it is not found, a local GF(2^128) implementation is
%      used. Keep encryption/decryption using the same offset convention.

    % ---------------------- input checks ----------------------
    if nargin < 7 || isempty(P) % nargin = number of input arguments
        P = zeros(1, n-k, 'uint8');
    end

    validateattributes(n, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(k, {'numeric'}, {'scalar','integer','positive','<',n});
    validateattributes(taglen, {'numeric'}, {'scalar','integer','positive','<=',n});

    if mod(n,8) ~= 0
        error('n must be a multiple of 8 because AES/Cipher is byte based.');
    end
    if length(P) ~= n-k
        error('P must contain exactly n-k bits.');
    end

    M = uint8(M(:).');
    P = uint8(P(:).');
    if any(M ~= 0 & M ~= 1) || any(P ~= 0 & P ~= 1)
        error('M and P must be binary bit vectors containing only 0 and 1.');
    end

    BLOCKLEN = n/8;
    ell_len = ceil(log2(n));        % 7 when n = 128
    if ell_len > k
        error('k must be at least ceil(log2(n)).');
    end

    % ---------------------- initialize L and H ----------------------
    % ------------------------------------------------------------
    % L <- E_K(0^n)
    % ------------------------------------------------------------
    zero_bytes = zeros(1, BLOCKLEN, 'uint8');
    
    zero_hex = upper(reshape(dec2hex(zero_bytes, 2).', 1, []));
    L_hex = Cipher(key_hex, zero_hex);
    L = uint8(hex2dec(reshape(L_hex, 2, []).')).';
    
    % ------------------------------------------------------------
    % H <- E_K(N xor L)
    % ------------------------------------------------------------
    N_bytes = nonce_to_bytes(N, n);
    
    H_input = bitxor(N_bytes, L);
    H_input_hex = upper(reshape(dec2hex(H_input, 2).', 1, []));
    H_hex = Cipher(key_hex, H_input_hex);
    H = uint8(hex2dec(reshape(H_hex, 2, []).')).';
    
    % ---------------------- partition M ----------------------
    mlen_bits = length(M);
    mbar = max(1, ceil(mlen_bits/k));

    % Regular blocks are M[1],...,M[mbar-1]. The last block is handled by
    % the special pseudocode branch below.
    regular_count = mbar - 1;
    Checksum = zeros(1, BLOCKLEN, 'uint8');

    % Determine if the final pseudocode branch creates one or two blocks.
    last_start = regular_count*k + 1;
    if last_start <= mlen_bits
        M_last = M(last_start:end);
    else
        M_last = uint8([]);
    end
    n_last = length(M_last);

    split_last = (n_last > k - ell_len);
    num_cipher_blocks = mbar + split_last;
    Z = rocb_offsets(L, H, num_cipher_blocks);

    C_out = zeros(1, num_cipher_blocks*n, 'uint8');
    C_bits_idx = 1;

    % ---------------------- regular blocks ----------------------
    for i = 1:regular_count
        block_bits = M((i-1)*k + 1 : i*k);
        B = [block_bits, P];

        [C_block_bits, B_bytes] = encrypt_rocb_block(key_hex, B, Z(i,:));
        C_out(C_bits_idx:C_bits_idx+n-1) = C_block_bits;
        C_bits_idx = C_bits_idx + n;

        Checksum = bitxor(Checksum, B_bytes);
    end

    % ---------------------- special last block(s) ----------------------
    P_tilde = [uint8(1 - P(1)), P(end-(n-k-1)+1:end)];
    L_x_inv = gf128_mul_x_inv(L);

    if ~split_last
        % n_last <= k - ell_len
        n_rand = k - ell_len - n_last;
        P_rand = uint8(randi([0,1], 1, n_rand));
        len_bits = uint_to_bits(n_rand, ell_len);

        B = [M_last, P_rand, len_bits, P_tilde];
        Z_last = bitxor(Z(mbar,:), L_x_inv);   % B xor (L*x^-1) xor Z[mbar]

        % [C_block_bits, B_bytes] = encrypt_rocb_block(key_hex, B, Z_last);
        [C_block_bits, B_bytes] =  encrypt_rocb_last_block(key_hex, B, Z_last, Z(mbar,:));
        C_out(C_bits_idx:C_bits_idx+n-1) = C_block_bits;
        Checksum = bitxor(Checksum, B_bytes);

        last_tag_offset_index = mbar;
        B_last = uint8([]);
    else
        % n_last > k - ell_len. Split M[mbar] into two encoded blocks.
        n_next = n_last - (k - ell_len);
        M_first = M_last(1:k-ell_len);
        M_second = M_last(k-ell_len+1:end);

        % First split block: no random padding, length field is 0^ell_len.
        zero_len_bits = zeros(1, ell_len, 'uint8');
        B = [M_first, zero_len_bits, P_tilde];

        [C_block_bits, B_bytes] = encrypt_rocb_block(key_hex, B, Z(mbar,:));
        C_out(C_bits_idx:C_bits_idx+n-1) = C_block_bits;
        C_bits_idx = C_bits_idx + n;
        Checksum = bitxor(Checksum, B_bytes);

        % Second split block: remaining bits + random padding + length.
        n_rand2 = k - ell_len - n_next;
        P_rand2 = uint8(randi([0,1], 1, n_rand2));
        len_bits2 = uint_to_bits(n_rand2, ell_len);
        B_last = [M_second, P_rand2, len_bits2, P_tilde];

        Z_last_masked = bitxor(Z(mbar+1,:), L_x_inv);
        % [C_block_bits, B_last_bytes] = encrypt_rocb_block(key_hex, B_last, Z_last_masked);
        [C_block_bits, B_last_bytes] = encrypt_rocb_last_block(key_hex, B_last, Z_last_masked, Z(mbar+1,:));
        C_out(C_bits_idx:C_bits_idx+n-1) = C_block_bits;
        Checksum = bitxor(Checksum, B_last_bytes);

        last_tag_offset_index = mbar + 1;
    end

    % ---------------------- tag ----------------------
    tag_input = bitxor(Checksum, Z(last_tag_offset_index,:));
    T_bytes = aes_enc_bytes(key_hex, tag_input);
    T_bits = bytes_to_bits(T_bytes);
    T_out = T_bits(1:taglen);

    % ---------------------- optional debug info ----------------------
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
        info.P_tilde = P_tilde;
        info.B_last = B_last;
        info.Checksum = Checksum;
        info.last_tag_offset_index = last_tag_offset_index;
    end
end

% =====================================================================
% Local helper functions
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

function [C_bits, B_bytes] = encrypt_rocb_last_block(key_hex, B_bits, Z_with_x_bytes, Z_bytes)
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

function out_bytes = aes_enc_bytes(key_hex, in_bytes)
    in_hex = upper(reshape(dec2hex(uint8(in_bytes), 2).', 1, []));
    out_hex = Cipher(key_hex, in_hex);
    out_bytes = uint8(hex2dec(reshape(out_hex, 2, []).')).';
end

function bytes = bits_to_bytes(bits)
    bits = uint8(bits(:).');
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
    % Prefer the existing project helper, if it exists, so encryption remains
    % compatible with the corresponding decryptor.
    if exist('compute_rocb_offsets', 'file') == 2
        Z_existing = compute_rocb_offsets(L, H, num_offsets);
        Z = uint8(Z_existing);
        return;
    end

    % Fallback: Z[i] = gamma_i*L xor H, with gamma_i represented by the
    % standard binary reflected Gray-code coefficient. If your thesis/code
    % defines gamma_i differently, replace only rocb_gamma_times_L().
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
    % Multiply a 128-bit field element by a small nonnegative integer whose
    % binary representation is interpreted as a GF(2) polynomial.
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
    % OCB/GCM-style doubling in GF(2^128) with polynomial
    % x^128 + x^7 + x^2 + x + 1. x is MSB-first byte order.
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
    % Inverse of gf128_mul_x(), i.e., multiplication by x^{-1}.
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
