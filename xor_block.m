function dst = xor_block(src1, src2)
    % XOR two 16-byte (128-bit) blocks
    dst = bitxor(uint8(src1), uint8(src2));
end