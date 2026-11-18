# R-OCB-MATLAB

MATLAB implementation of R-OCB, including encryption/decryption, GRAND/ORBGRAND decoding, noisy-channel simulations, CRC variants, baselines, and plotting scripts used in the final paper.

## R-OCB: Reliable Offset Codebook

This repository contains the MATLAB implementation and simulation code used in the final paper:

**R-OCB: Reliable Offset Codebook**  
Noga Cohen, Technion – Israel Institute of Technology

The implementation evaluates R-OCB, a block-cipher-based mode that combines encryption with structured redundancy for reliable communication over noisy channels.

## GRAND / ORBGRAND

The GRAND and ORBGRAND implementations and several comparison-code constructions are based on the public GRAND-MATLAB repository:

https://github.com/kenrduffy/GRAND-MATLAB

The R-OCB implementation and simulation code corresponding to the results reported in the final paper are archived in this repository. The exact matrix file used for the retained CA-Polar (n=128), (k=116) baseline is identified in the final paper.
