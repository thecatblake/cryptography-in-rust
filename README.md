# Cryptography Implementation Roadmap

An ongoing personal journey to deepen my understanding of modern cryptography by implementing everything from the mathematical foundations to advanced zero-knowledge proof systems, Ethereum internals, and secure multi-party computation — from scratch.

**Goals:** contribute to open-source cryptography projects, read research papers fluently, and implement new cryptographic protocols.

---

# My View on AI (for Learners)
The more you rely on AI, the less often you discover what you don't know—your unknown unknowns.
The less you rely on AI, the slower your work becomes, but the more opportunities you have to uncover gaps in your understanding.

# Phase 0 — Core Mathematics

## BigInt

* [x] Addition
* [x] Subtraction
* [x] Multiplication
* [x] Division
* [x] Modulo
* [x] Comparison
* [x] GCD / Extended GCD
* [x] Modular Inverse
* [x] Montgomery Reduction
* [ ] Barrett Reduction
* [ ] Modular Exponentiation
* [ ] Miller–Rabin Primality Test
* [ ] Constant-time comparison and conditional select
* [ ] Karatsuba multiplication
* [ ] ◇ Toom–Cook multiplication
* [ ] Fixed-width limb representation (`u64` limbs, carry handling)

---

# Phase 1 — Finite Fields

## Prime Field (`Fp`)

* [x] Addition
* [x] Subtraction
* [x] Multiplication
* [x] Division
* [x] Negation
* [x] Inversion
* [x] Exponentiation (`pow`)
* [x] Squaring
* [ ] Square Root (Tonelli–Shanks)
* [ ] Legendre symbol / quadratic residuosity
* [ ] Montgomery form arithmetic (REDC)
* [ ] Batch Inversion
* [ ] Random field element sampling
* [ ] Serialization (canonical byte encoding)

## Small Fields ⚑

*Modern STARKs live here, not in 256-bit fields. Cheaper to implement than `Fp` and they make Phase 5 tractable.*

* [x] Goldilocks (`p = 2^64 − 2^32 + 1`)
* [x] BabyBear (`p = 2^31 − 2^27 + 1`)
* [x] Mersenne31 (`p = 2^31 − 1`)
* [ ] Binary fields / tower field construction (`GF(2^k)`)
* [ ] Field extension for small fields (soundness needs a large enough extension)

## Extension Fields

* [x] `Fp2`
* [x] `Fp6`
* [x] `Fp12`
* [x] Frobenius endomorphism (generic on `QuadExt`/`CubicExt`; `Fp2`/`Fp6` coefficients computed at compile time, `Fp12`'s verified offline)
* [ ] Towering strategy and multiplication cost analysis

## Field Infrastructure

* [x] Generic `Field` / `PrimeField` trait design (`Field` trait + `QuadExt`/`CubicExt` generic extension towers; `Fp2`/`Fp6` are now derived instantiations)
* [ ] FFT-friendly field (two-adicity, root-of-unity generation)
* [ ] Multiplicative subgroup / coset construction

---

# Phase 2 — Elliptic Curves

## Fundamentals

* [x] Curve Definition (short Weierstrass)
* [x] Point Representation (affine)
* [x] Identity Point
* [x] Point Addition
* [x] Point Doubling
* [x] Scalar Multiplication (double-and-add)
* [x] Point Validation (on-curve)
* [ ] Point Validation (subgroup check)
* [ ] Serialization / Deserialization (compressed + uncompressed)

## Coordinate Systems

* [ ] Projective coordinates
* [x] Jacobian coordinates (short Weierstrass: inversion-free add/double, `to_affine`/`from_affine`, on-curve validation)
* [x] Extended twisted Edwards coordinates (unified add + dedicated double, `to_affine`/`from_affine`, on-curve validation)
* [ ] Cost comparison table (M/S/A counts per operation)

## Curve Forms

* [x] Montgomery form (curve definition + unified chord-and-tangent addition and double-and-add scalar multiplication; Curve25519 itself not wired up yet)
* [x] Twisted Edwards curve definition
* [x] Twisted Edwards point arithmetic (unified addition + double-and-add scalar multiplication; Ed25519 itself not wired up yet)
* [ ] In-circuit friendly curves (Baby Jubjub, Jubjub, Grumpkin)

## Advanced Scalar Multiplication

* [x] Montgomery ladder (constant-time and variable-time, plus the original double-and-add, selected at compile time via `elliptic-curve`'s mutually-exclusive `scalar-mul-double-and-add` / `scalar-mul-ladder-variable` / `scalar-mul-ladder-constant` Cargo features)
* [ ] Windowed scalar multiplication
* [ ] wNAF representation
* [ ] GLV endomorphism decomposition (secp256k1, BLS12-381)
* [ ] Multi-Scalar Multiplication (Pippenger / bucket method) ⚑
* [ ] Precomputed table generation

## Hashing to Curves

* [ ] Hash-to-field (RFC 9380)
* [ ] Simplified SWU map
* [ ] Elligator 2
* [ ] Clearing the cofactor
* [ ] Full `hash_to_curve` for BLS12-381 G1/G2

## Named Curves

* [ ] secp256k1
* [ ] P-256
* [ ] Curve25519 / Ed25519
* [ ] BN254
* [ ] BLS12-381
* [ ] Pasta curves (Pallas / Vesta) — 2-cycle
* [ ] ◇ BLS12-377 / BW6-761

---

# Phase 3 — Symmetric Primitives & Hashing

## Traditional Hashes

* [ ] SHA-2 (SHA-256, SHA-512)
* [ ] SHA-3
* [ ] Keccak-f[1600] and Keccak-256 (Ethereum variant)
* [ ] BLAKE2 / BLAKE3
* [ ] RIPEMD-160
* [ ] Merkle–Damgård vs. sponge construction notes

## Keyed / Derived

* [ ] HMAC
* [ ] HKDF
* [ ] PRF abstraction
* [ ] ◇ Argon2 / scrypt

## Symmetric Encryption

*Needed for MPC transport, hybrid encryption, and keystore formats.*

* [ ] AES block cipher
* [ ] AES-GCM
* [ ] ChaCha20
* [ ] Poly1305
* [ ] ChaCha20-Poly1305 AEAD

## Algebraic (ZK-Friendly) Hashes

* [ ] Poseidon
* [ ] Poseidon2
* [ ] Rescue / Rescue-Prime
* [ ] Pedersen Hash
* [ ] ◇ Griffin
* [ ] ◇ Monolith
* [ ] ◇ Skyscraper
* [ ] Sponge / duplex construction over a field

## Merkle Trees

* [ ] Binary Merkle Tree
* [ ] Merkle Proof generation and verification
* [ ] Sparse Merkle Tree
* [ ] Incremental Merkle Tree
* [ ] Multi-proofs / batch openings
* [ ] Domain separation and second-preimage resistance
* [ ] Merkle Mountain Range

---

# Phase 4 — Classical Public-Key Cryptography

## Randomness

* [ ] CSPRNG interface (`rand_core`-style traits)
* [ ] Rejection sampling for uniform field elements
* [ ] Deterministic nonce generation (RFC 6979)
* [ ] Hash-to-field (RFC 9380)

## Signatures

* [ ] Schnorr Signature
* [ ] EdDSA (Ed25519, including RFC 8032 edge cases)
* [ ] ECDSA
* [ ] ECDSA with public key recovery (`ecrecover`), low-`s` normalization
* [ ] BLS Signature (sign, verify, aggregate)
* [ ] Proof of Possession / rogue-key attack mitigation
* [ ] Batch verification

## Encryption & Key Exchange

* [ ] Diffie–Hellman
* [ ] ECDH / X25519
* [ ] RSA (keygen, encrypt, sign, PSS/OAEP padding)
* [ ] ElGamal
* [ ] Paillier (additively homomorphic)
* [ ] ◇ ECIES

## Secret Sharing

* [ ] Shamir's Secret Sharing
* [ ] Additive secret sharing
* [ ] Feldman VSS
* [ ] Pedersen VSS

## Commitment Schemes

* [ ] Pedersen Commitment
* [ ] Vector Pedersen Commitment
* [ ] Hash-based commitment
* [ ] Binding / hiding property analysis

---

# Phase 5 — Polynomials, Pairings & Proof Systems

## Polynomials

* [ ] Polynomial Representation (coefficient form)
* [ ] Polynomial Evaluation (Horner)
* [ ] Polynomial arithmetic (add, mul, div, mod)
* [ ] Polynomial Interpolation
* [ ] Lagrange Basis
* [ ] Vanishing polynomial over a domain
* [ ] FFT / NTT ⚑
* [ ] Inverse FFT (IFFT)
* [ ] Coset FFT
* [ ] Multilinear extensions and the `eq` polynomial ⚑
* [ ] Univariate ↔ multilinear conversions

## Pairings

* [ ] Miller Loop
* [ ] Final Exponentiation
* [ ] BN254 Pairing
* [ ] BLS12-381 Pairing
* [ ] Multi-pairing / batch pairing check
* [ ] ◇ Optimal ate pairing optimizations

## Interactive Proof Machinery

* [ ] Sumcheck protocol ⚑ *(do this before Groth16 — it's simpler and unlocks GKR, HyperNova, Lasso)*
* [ ] Fiat–Shamir Transform
* [ ] Transcript (with domain separation)
* [ ] Sigma protocols / proofs of knowledge
* [ ] OR-composition
* [ ] GKR
* [ ] Soundness accounting and security-bit budgeting
* [ ] Grinding / proof-of-work in transcripts

## Arithmetization

* [ ] Rank-1 Constraint System (R1CS)
* [ ] Quadratic Arithmetic Program (QAP)
* [ ] Constraint System abstraction
* [ ] Witness generation
* [ ] Minimal circuit DSL
* [ ] AIR (Algebraic Intermediate Representation)
* [ ] RAP (AIR with challenges / preprocessed columns)
* [ ] Plonkish: custom gates, copy constraints, grand-product permutation argument
* [ ] ◇ CCS (generalizes R1CS / Plonkish / AIR)

## Lookup Arguments

* [ ] Plookup
* [ ] LogUp
* [ ] ◇ LogUp-GKR
* [ ] ◇ Lasso

## Circuit Gadgets

* [ ] Range checks
* [ ] Bit decomposition
* [ ] Non-native (foreign) field arithmetic
* [ ] In-circuit elliptic curve operations
* [ ] In-circuit hashing (Poseidon)
* [ ] In-circuit Merkle proof verification

## Polynomial Commitment Schemes

* [ ] KZG Commitment
* [ ] Batch / multi-point KZG openings
* [ ] Inner Product Argument (IPA) Commitment
* [ ] Bulletproofs
* [ ] FRI ⚑
* [ ] DEEP-FRI, batched FRI
* [ ] ◇ Ligero / Brakedown
* [ ] ◇ Basefold
* [ ] ◇ WHIR
* [ ] ◇ Binius (binary-field, tower construction)
* [ ] ◇ Zeromorph / Mercury / Dory (multilinear from KZG)

## Proof Systems

* [ ] STARK over Goldilocks *(suggested first system — no pairings, no trusted setup)*
* [ ] Groth16
* [ ] PLONK
* [ ] ◇ Halo2 *(highest complexity-to-insight ratio — skip unless specifically needed)*
* [ ] ◇ Circle STARKs over Mersenne31
* [ ] Groth16 verifier in Solidity (gas-optimized)

## Folding & Recursion

* [ ] Accumulation schemes (abstract)
* [ ] Nova
* [ ] ◇ SuperNova
* [ ] ◇ HyperNova
* [ ] ◇ ProtoStar / ProtoGalaxy
* [ ] Two-cycle curves in practice (Pasta, BN254/Grumpkin)
* [ ] ◇ CycleFold
* [ ] In-circuit verifier for your own SNARK ⚑ *(the real test)*
* [ ] Proof aggregation and batch verification

## zkVM

* [ ] RISC-V instruction encoding and decoding
* [ ] Memory checking (offline / permutation-based)
* [ ] Register and program-counter constraints
* [ ] Continuations / chunked proving
* [ ] ◇ Precompiles for hashing and field ops

---

# Phase 6 — Ethereum

## Encoding & State

* [ ] RLP encode / decode
* [ ] SSZ encode / decode, Merkleization, `hash_tree_root`
* [ ] Merkle Patricia Trie (branch/extension/leaf nodes, hex-prefix encoding)
* [ ] MPT proofs (inclusion and exclusion)
* [ ] Storage / transaction / receipt tries
* [ ] ABI encoding / decoding (static, dynamic, packed)
* [ ] EIP-712 typed structured data hashing
* [ ] ◇ Verkle trees + Banderwagon / IPA
* [ ] ◇ Binary trie with Poseidon (the alternative proposal)

## Accounts & Transactions

* [ ] Address derivation from public key
* [ ] CREATE / CREATE2 address computation
* [ ] EIP-155 chain-id replay protection
* [ ] Legacy transactions
* [ ] EIP-2930 (access lists)
* [ ] EIP-1559 (fee market, base fee calculation)
* [ ] EIP-4844 (blob transactions)
* [ ] EIP-7702 (set EOA code)
* [ ] ◇ ERC-4337 UserOperation validation
* [ ] BIP-32 / BIP-39 / BIP-44 key derivation
* [ ] Keystore (Web3 Secret Storage) encrypt / decrypt

## Execution Layer

* [ ] EVM interpreter: stack, memory, calldata, storage
* [ ] Full opcode set
* [ ] Gas metering (including EIP-2929 warm/cold access)
* [ ] Call frames: CALL, DELEGATECALL, STATICCALL
* [ ] Revert semantics and journaling / rollback
* [ ] Precompiles: `ecrecover`, `sha256`, `ripemd160`, identity, `modexp`
* [ ] Precompiles: BN254 add / mul / pairing
* [ ] Precompiles: `blake2f`
* [ ] Precompile: KZG point evaluation (EIP-4844)
* [ ] State transition function
* [ ] Block and header validation
* [ ] Withdrawal processing

## Consensus Layer

* [ ] BLS12-381 signature aggregation (min-pubkey-size scheme)
* [ ] Beacon state and block processing
* [ ] Attestations
* [ ] Sync committees
* [ ] Committee shuffling (`compute_shuffled_index`)
* [ ] LMD-GHOST fork choice
* [ ] Casper FFG finality
* [ ] Light client: sync committee proof verification
* [ ] 4844 blobs: KZG commitment, versioned hash, blob proof verification
* [ ] ◇ Data availability sampling / PeerDAS (EIP-7594)

## Networking ◇

*Optional, but where a lot of real-world bugs live.*

* [ ] devp2p RLPx handshake
* [ ] discv5
* [ ] Snap sync / range proofs
* [ ] libp2p gossipsub (consensus layer)

---

# Phase 7 — Multi-Party Computation

## Foundations

* [ ] Security model taxonomy: semi-honest vs. malicious, honest vs. dishonest majority
* [ ] Simulation-based security definitions — read the proofs, not just the protocols ⚑
* [ ] Ideal/real paradigm and functionality descriptions
* [ ] Composition (UC) basics

## Secret Sharing (MPC-specific)

* [ ] Replicated secret sharing
* [ ] Packed / ramp secret sharing
* [ ] Distributed Key Generation (Pedersen DKG)
* [ ] New-DKG / secure variants
* [ ] Proactive refresh and resharing

## Oblivious Transfer

* [ ] Base OT (Naor–Pinkas)
* [ ] Base OT (Chou–Orlandi "simplest OT")
* [ ] OT extension: IKNP (semi-honest)
* [ ] OT extension: KOS (malicious)
* [ ] ◇ SoftSpokenOT
* [ ] Correlated OT / random OT
* [ ] VOLE
* [ ] ◇ Silent OT (Ferret), pseudorandom correlation generators
* [ ] Function secret sharing / Distributed Point Functions

## Circuit-Based MPC

* [ ] Yao's garbled circuits (textbook)
* [ ] Point-and-permute
* [ ] Free-XOR
* [ ] Half-gates
* [ ] ◇ Three-halves garbling
* [ ] GMW (boolean)
* [ ] GMW (arithmetic)
* [ ] ◇ BMR (multiparty garbling)
* [ ] Beaver triples and the preprocessing model
* [ ] BGW / Shamir-based multiplication with degree reduction
* [ ] SPDZ: authenticated shares, MAC check, sacrifice
* [ ] ◇ MASCOT / Overdrive triple generation
* [ ] Cut-and-choose for malicious security

## Threshold Cryptography

* [ ] Threshold Schnorr — FROST
* [ ] ◇ ROAST
* [ ] Threshold BLS
* [ ] Threshold ECDSA — Lindell
* [ ] ◇ Threshold ECDSA — GG18 / GG20
* [ ] ◇ Threshold ECDSA — CGGMP21
* [ ] Threshold decryption (ElGamal, Paillier)

## Applications

* [ ] Private Set Intersection
* [ ] Oblivious PRF
* [ ] ◇ OPAQUE (password-authenticated key exchange)
* [ ] Trusted-setup MPC ceremony: contribution + verification (powers of tau)
* [ ] ◇ Private information retrieval

---

# Phase 8 — Lattices & Post-Quantum ◇

*Parallel track. Kyber and Dilithium are more approachable than Halo2 — precise specs, exhaustive test vectors, and NTT is close to work already done in Phase 5.*

## Foundations

* [ ] LWE / RLWE / Module-LWE
* [ ] NTT and negacyclic convolution
* [ ] Gaussian sampling
* [ ] Centered binomial distribution sampling
* [ ] Lattice reduction basics (LLL) for parameter intuition

## Standardized Schemes

* [ ] ML-KEM (Kyber)
* [ ] ML-DSA (Dilithium)
* [ ] ◇ Falcon
* [ ] SLH-DSA (SPHINCS+)
* [ ] XMSS / hash-based signatures

## Fully Homomorphic Encryption

* [ ] BFV
* [ ] BGV
* [ ] CKKS (approximate arithmetic)
* [ ] Relinearization and key switching
* [ ] Modulus switching and noise budget analysis
* [ ] ◇ TFHE and programmable bootstrapping
* [ ] ◇ Threshold FHE decryption

---

# Cross-Cutting: Testing & Correctness

* [ ] Unit tests with known-answer vectors
* [ ] Property Testing (`proptest` / `quickcheck`)
* [ ] Fuzz Testing (`cargo-fuzz`)
* [ ] Differential fuzzing against reference implementations
* [ ] Test vectors from RFC 9380, RFC 6979, RFC 8032, NIST CAVP
* [ ] Ethereum `execution-spec-tests` as fixtures
* [ ] Ethereum `consensus-spec-tests` as fixtures
* [ ] Cross-check against `arkworks`, `blst`, `revm`, `geth`
* [ ] Constant-Time Verification (`subtle`, branch-free select) -- the `scalar-mul-ladder-constant` feature gets the ladder's *shape* right (fixed iteration count, arithmetic `cswap` instead of branching on the coordinates) but still branches on the bit when turning it into a field element (see `elliptic-curve/src/ladder.rs`); closing that gap needs a `Field`-level branch-free select primitive, which is what this item is for
* [ ] `dudect` / statistical timing analysis
* [ ] ◇ Formal verification of one field or curve routine (hacspec, Kani, Creusot)
* [ ] Negative tests: malformed inputs, invalid points, subgroup attacks
* [ ] Soundness tests: proofs of false statements must fail

---

# Cross-Cutting: Performance

* [x] Benchmarking harness (`criterion`)
* [x] Montgomery field arithmetic
* [ ] SIMD optimization (AVX2 / AVX-512 / NEON)
* [ ] Parallel MSM
* [ ] Parallel FFT
* [ ] Parallel Merkle tree construction
* [ ] ◇ Assembly optimization for field multiplication
* [ ] ◇ GPU acceleration (MSM / NTT)
* [ ] Memory layout and cache-friendly representations
* [ ] `no_std` support
* [ ] WASM target support
* [ ] Profiling (`perf`, `flamegraph`)

### Benchmark Results

Run via `cargo bench --workspace`. Median timings from `criterion`; local-machine numbers, not portable across hardware — useful for relative comparison between backends, not absolute performance claims.

| op | BabyBear (Montgomery, u32) | Goldilocks (native, u64) | Mersenne31 (native, u32) | `Fp<U256>` DefaultBackend | `Fp<U256>` MontBackend |
|---|---|---|---|---|---|
| add | 629 ps | 552 ps | 631 ps | 1.60 ns | 1.60 ns |
| sub | 563 ps | 556 ps | 561 ps | 1.84 ns | 1.77 ns |
| mul | 737 ps | 660 ps | 625 ps | 632 ns | 20.4 ns |
| square | 527 ps | 592 ps | 435 ps | 601 ns | 19.9 ns |
| neg | 336 ps | 341 ps | 342 ps | 1.60 ns | 1.60 ns |
| inverse | 88.5 ns | 229 ns | 76.0 ns | 877 ns | 921 ns |
| pow | 81.6 ns | 98.9 ns | 69.8 ns | 408 µs | 8.64 µs |

Montgomery reduction cuts `Fp<U256>` `mul` from 632ns to 20.4ns (~31x) by replacing a division with multiply+shift. BabyBear and Goldilocks are a further order of magnitude faster than `Fp<U256>` since they operate on a single machine word instead of a multi-limb bigint.

### Extension Field Benchmark Results (Fp2 / Fp6 / Fp12)

Run via `cargo bench --workspace`. Median timings from `criterion`; `Fp<U256>` uses `MontBackend` over the secp256k1 base field prime (see `field_ext_ops.rs`); same local-machine caveat as above.

**Fp2**

| op | BabyBear | Goldilocks | Mersenne31 | `Fp<U256>` MontBackend |
|---|---|---|---|---|
| add | 843 ps | 803 ps | 846 ps | 3.19 ns |
| sub | 841 ps | 795 ps | 837 ps | 5.92 ns |
| mul | 2.48 ns | 3.54 ns | 2.53 ns | 99.3 ns |
| square | 2.37 ns | 3.30 ns | 2.39 ns | 95.9 ns |
| neg | 585 ps | 567 ps | 560 ps | 3.43 ns |
| inverse | 112 ns | 249 ns | 89.6 ns | 1.04 µs |

**Fp6**

| op | BabyBear | Goldilocks | Mersenne31 | `Fp<U256>` MontBackend |
|---|---|---|---|---|
| add | 2.23 ns | 2.82 ns | 2.12 ns | 16.0 ns |
| sub | 2.24 ns | 2.68 ns | 2.16 ns | 19.8 ns |
| mul | 26.5 ns | 34.8 ns | 29.2 ns | 807 ns |
| square | 26.4 ns | 34.4 ns | 29.2 ns | 806 ns |
| neg | 1.51 ns | 1.69 ns | 1.43 ns | 11.0 ns |
| inverse | 174 ns | 310 ns | 167 ns | 2.06 µs |

**Fp12**

| op | BabyBear | Goldilocks | Mersenne31 | `Fp<U256>` MontBackend |
|---|---|---|---|---|
| add | 4.34 ns | 6.82 ns | 4.38 ns | 45.9 ns |
| sub | 4.33 ns | 4.76 ns | 4.21 ns | 58.1 ns |
| mul | 119 ns | 153 ns | 170 ns | 3.30 µs |
| square | 118 ns | 152 ns | 168 ns | 3.24 µs |
| neg | 2.84 ns | 3.26 ns | 2.91 ns | 31.7 ns |
| inverse | 320 ns | 491 ns | 388 ns | 6.37 µs |

`mul`/`square`/`inverse` grow roughly 3x per tower level (Fp2 → Fp6 → Fp12) on the small fields, tracking the Karatsuba mul cost (3 base muls at Fp2, 6 at Fp6, 3 Fp6-muls at Fp12) plus the extra field inversion each `inverse` call chains through. `Fp<U256>`'s multi-limb `mul` dominates every tower level, same as it does for the base field above.

### Elliptic Curve Benchmark Results

Run via `cargo bench -p elliptic-curve` (`elliptic_curve_ops.rs`). Median timings from `criterion`; all three curve forms run over the same `Fp<U256>` `MontBackend` field (secp256k1's base prime, reused only for a realistic 256-bit modulus — the curve constants and points are arbitrary, not secp256k1's own, since point `add`/`Mul` never call `validate()`), the backend real curve code would actually use for its cheaper `mul` (see the `Fp<U256>` rows in the field benchmark table above); `scalar_mul` uses a full-width 256-bit scalar. Same local-machine caveat as above.

| op | Short Weierstrass | Twisted Edwards | Montgomery |
|---|---|---|---|
| add | 707 ns | 3.04 µs | 738 ns |
| double | 1.49 µs | — | 1.56 µs |
| scalar_mul | 821 µs | 1.79 ms | 858 µs |

Twisted Edwards has no separate `double`: its addition law is unified (see `twisted_edwards.rs`), so doubling reuses the `add` code path rather than a cheaper tangent-line formula — that's also why its `add` costs roughly 2 inversions' (via 2 divisions) worth more than short Weierstrass/Montgomery's chord-and-tangent `add`, which only needs 1. `scalar_mul` (double-and-add over a 256-bit scalar) tracks each form's `add`/`double` cost roughly linearly, since it's ~256 doublings plus up to 256 adds. Switching the field backend from `DefaultBackend` to `MontBackend` (measured directly: `add` 2.28 µs → 707 ns, `double` 3.73 µs → 1.49 µs, `scalar_mul` 2.36 ms → 821 µs on Short Weierstrass) cuts every curve op by roughly 2.5–3x, in line with `MontBackend`'s cheaper `mul` in the base-field table above dominating the handful of field multiplications each curve op chains through.

---

# Cross-Cutting: Security Hygiene

* [ ] No secret-dependent branches
* [ ] No secret-dependent memory indexing
* [ ] Zeroization of secret material (`zeroize`)
* [ ] Domain separation everywhere (hashes, transcripts, signatures)
* [ ] Nonce misuse resistance review
* [ ] Subgroup and small-order point checks on all deserialization
* [ ] Malleability review (signatures, encodings)
* [ ] Documented security assumptions per module
* [ ] Threat model written down before implementing each protocol

---
