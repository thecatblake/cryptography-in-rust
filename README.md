For contributing to open-source cryptography projects, understanding research papers, and implementing new cryptographic protocols, I created this roadmap.

This project is an ongoing personal journey to deepen my understanding of modern cryptography by implementing everything from the mathematical foundations to advanced zero-knowledge proof systems from scratch.

## Core Mathematics

* [ ] Implement `BigInt`

  * [x] Addition
  * [x] Subtraction
  * [x] Multiplication
  * [x] Division
  * [x] Modulo
  * [x] Comparison
  * [x] GCD / Extended GCD
  * [ ] Modular Inverse
  * [ ] Montgomery Reduction
  * [ ] Barrett Reduction
  * [ ] Modular Exponentiation
  * [ ] Miller–Rabin Primality Test

## Finite Fields (`field`)

* [ ] Prime Field (`Fp`)

  * [ ] Addition
  * [ ] Subtraction
  * [ ] Multiplication
  * [ ] Division
  * [ ] Negation
  * [ ] Inversion
  * [ ] Exponentiation (`pow`)
  * [ ] Squaring
  * [ ] Square Root (Tonelli–Shanks)

* [ ] Extension Fields

  * [ ] `Fp2`
  * [ ] `Fp6`
  * [ ] `Fp12`

* [ ] FFT-Friendly Field

* [ ] Batch Inversion

## Elliptic Curves (`ecc`)

* [ ] Curve Definition
* [ ] Point Representation
* [ ] Identity Point
* [ ] Point Addition
* [ ] Point Doubling
* [ ] Scalar Multiplication
* [ ] Windowed Scalar Multiplication
* [ ] Multi-Scalar Multiplication (MSM)
* [ ] Serialization / Deserialization
* [ ] Point Validation
* [ ] Hash-to-Curve

## Hash Functions

* [ ] SHA-2
* [ ] SHA-3
* [ ] Keccak
* [ ] Poseidon
* [ ] Rescue
* [ ] Pedersen Hash

## Digital Signatures

* [ ] Schnorr Signature
* [ ] EdDSA
* [ ] ECDSA
* [ ] BLS Signature

## Commitment Schemes

* [ ] Pedersen Commitment
* [ ] KZG Commitment
* [ ] Inner Product Argument (IPA) Commitment

## Merkle Trees

* [ ] Binary Merkle Tree
* [ ] Sparse Merkle Tree
* [ ] Merkle Proof
* [ ] Incremental Merkle Tree

## Pairings

* [ ] Miller Loop
* [ ] Final Exponentiation
* [ ] BN254 Pairing
* [ ] BLS12-381 Pairing

## Polynomials

* [ ] Polynomial Representation
* [ ] Polynomial Evaluation
* [ ] Polynomial Interpolation
* [ ] Lagrange Basis
* [ ] FFT
* [ ] Inverse FFT (IFFT)
* [ ] Polynomial Commitments

## SNARK Foundations

* [ ] Rank-1 Constraint System (R1CS)
* [ ] Quadratic Arithmetic Program (QAP)
* [ ] Fiat–Shamir Transform
* [ ] Transcript
* [ ] Constraint System

## Zero-Knowledge Proof Systems

* [ ] Groth16
* [ ] PLONK
* [ ] Halo2
* [ ] Nova
* [ ] STARK (FRI)

## Cryptographic Protocols

* [ ] Diffie–Hellman
* [ ] ElGamal
* [ ] Paillier
* [ ] RSA
* [ ] Threshold Signatures
* [ ] Shamir's Secret Sharing

## Testing

* [ ] Property Testing
* [ ] Fuzz Testing
* [ ] Benchmarking (`criterion`)
* [ ] Constant-Time Verification
* [ ] Test Vectors Against Existing Libraries

## Optimization

* [ ] SIMD Optimization
* [ ] Parallel MSM
* [ ] Parallel FFT
* [ ] Montgomery Field Arithmetic
* [ ] Assembly Optimization
* [ ] `no_std` Support
