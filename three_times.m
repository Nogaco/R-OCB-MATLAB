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
 
 function dst = three_times(src)
    % Multiply 128-bit block by 3 in GF(2^128): 3*X = 2*X ^ X
    t = two_times(src);
    dst = xor_block(t, src);
end
