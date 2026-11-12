%/*************************************************************************** 
% * two_times, three_times
% *
% * These routines manipulate the offsets which are used for pre- and
% * post-whitening of blockcipher invocations. The offsets represent
% * polynomials, and these routines multiply the polynomials by other
% * constant polynomials. Note that as an optimization two consecutive
% * invocations of "three_times" can be efficiently replaced:
% *                     3(3(X)) == (2(2(X))) xor X
% *
% ***************************************************************************/
 
 function dst = two_times(src)
    % Multiply 128-bit block by 2 in GF(2^128)
    src = uint8(src);
    dst = zeros(1, 16, 'uint8');
    carry = bitshift(src(1), -7);  % Get MSB of first byte

    for i = 1:15
        dst(i) = bitshift(src(i), 1) + bitget(src(i+1), 8);
    end

    dst(16) = bitxor(bitshift(src(16), 1), uint8(carry * hex2dec('87')));
end