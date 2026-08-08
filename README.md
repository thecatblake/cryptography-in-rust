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
* [ ] `Fp6`
* [ ] `Fp12`
* [ ] Frobenius endomorphism
* [ ] Towering strategy and multiplication cost analysis

## Field Infrastructure

* [ ] Generic `Field` / `PrimeField` trait design
* [ ] FFT-friendly field (two-adicity, root-of-unity generation)
* [ ] Multiplicative subgroup / coset construction

---

# Phase 2 — Elliptic Curves

## Fundamentals

* [ ] Curve Definition (short Weierstrass)
* [ ] Point Representation (affine)
* [ ] Identity Point
* [ ] Point Addition
* [ ] Point Doubling
* [ ] Scalar Multiplication (double-and-add)
* [ ] Point Validation (on-curve, subgroup check)
* [ ] Serialization / Deserialization (compressed + uncompressed)

## Coordinate Systems

* [ ] Projective coordinates
* [ ] Jacobian coordinates
* [ ] Extended twisted Edwards coordinates
* [ ] Cost comparison table (M/S/A counts per operation)

## Curve Forms

* [ ] Montgomery form (Curve25519)
* [ ] Twisted Edwards form (Ed25519)
* [ ] In-circuit friendly curves (Baby Jubjub, Jubjub, Grumpkin)

## Advanced Scalar Multiplication

* [ ] Montgomery ladder (constant-time)
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
* [ ] Constant-Time Verification (`subtle`, branch-free select)
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
