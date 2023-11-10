# TREnc-PPAT

This repository contains the companion code for the paper "An encryption mechanism for receipt-free and perfectly private verifiable elections" by Thi Van Thao Doan, Olivier Pereira, and Thomas Peters.

In this work, we present two efficient encryption mechanisms for verifiable voting protocols, adapted to the two most crucial tallying techniques:
- one mechanism is additively homomorphic and offers a polylogarithmic message space, suitable for homomorphic ballot aggregation;
- one mechanism is randomizable and supports the encryption of full-size group elements, making it compatible with a mixnet-based tallying process.

To assess the practicality of these mechanisms, we benchmark prototype implementations
of the key generation, encryption and verification algorithms. We rely on the
MIRACL Core Cryptographic Library [MIRACL](https://github.com/miracl/core) for the elliptic curve and pairing operations on the BN462 curve. 

## Source

The repository contains the following files:
- `tools.h`: contains related-object structures and shared functions between two constructions.
- `TREnc_PPAT_simple_ballots.c`: contains functions of key generation, encryption, verification, and decryption for our additively homomorphic mechanism (simple ballots).
- `TREnc_PPAT_complex_ballots.c`: contains functions of key generation, encryption, verification, and decryption for our randomizable mechanism (complex ballots).


## Benchmarks

The repository contains the following benchmarks:
- `benchmark_complexBallot.c`
- `benchmark_simpleBallot.c`

Each one of these benchmarks can be compiled using `gcc` command (tested with version 11.4.0).
```
$ gcc -O3 filename.c library/core.a -o filename
$ ./filename
```
**Output examples:**
```shell
Testing/Timing the scheme for complex ballots
KEYGEN RUNNING TIME             -       10 iterations      11.41 ms per iteration
ENCRYPTION RUNNING TIME         -       10 iterations     142.35 ms per iteration
VERIFICATION RUNNING TIME       -       10 iterations     537.83 ms per iteration
DECRYPTION IS CORRECT!
```

```shell
Testing/Timing the scheme for simple ballots
KEYGEN RUNNING TIME             -       10 iterations      14.24 ms per iteration
ENCRYPTION RUNNING TIME         -       10 iterations     154.99 ms per iteration
VERIFICATION RUNNING TIME       -       10 iterations     550.98 ms per iteration
DECRYPTION IS CORRECT!
```
## Acknowledgements
 The research depicted in this paper is funded by CYBEREXCELLENCE, a cyber security excellence project as part of the Walloon Region plan (CyberWal) in Belgium.
