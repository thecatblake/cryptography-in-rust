# Ethereum Roadmap Research (as of 2026-08)

A summary of Ethereum's current state and near-term plans, compiled from the Ethereum Foundation's official site (ethereum.org/roadmap), the EIP repository, Forkcast, Vitalik Buterin's public writing, and The Block. Written to be cross-referenced against `README.md`'s Phase 6 (Ethereum) and Phase 8 (Lattices & Post-Quantum).

---

## 1. The big picture: The Merge / Surge / Scourge / Verge / Purge / Splurge

Since the 2022 Merge, the roadmap has been narrated through five informal phases coined by Vitalik. Through 2025-2026, the actual forks (Pectra, Fusaka, Glamsterdam, Hegota...) implement these themes in overlapping fashion rather than one phase at a time.

| Phase | Goal | Current status (fork / EIP) |
|---|---|---|
| **The Merge** | PoW to PoS transition | Done (Paris, 2022-09-15) |
| **The Surge** | L2/blob scaling (data availability sampling, danksharding) | Proto-Danksharding shipped (Dencun, EIP-4844); PeerDAS shipped (Fusaka, EIP-7594); Full Danksharding still ahead |
| **The Scourge** | Removing MEV-driven centralization and censorship risk | ePBS planned (EIP-7732, Glamsterdam); FOCIL planned (EIP-7805, Hegota) |
| **The Verge** | Efficient verification / statelessness (originally Verkle trees) | Direction shifting: Verkle trees are being displaced by a **binary state tree (EIP-7864)** as the leading candidate |
| **The Purge** | Removing old history and unneeded state, protocol simplification | Partial history expiry shipped (2025-07); EIP-4444 (full history) and state expiry (research) still ahead |
| **The Splurge** | Wrapping up the above, general polish (gas repricing, UX, intrinsic gas cuts, etc.) | Most of Glamsterdam's smaller EIPs fall under this heading |

As of February 2026, Vitalik has added two more items to this list of "changes Ethereum needs to make while in flight": Lean Consensus (post-quantum consensus) and a VM overhaul (replacing the EVM with RISC-V or WASM).

**Sources**: ethereum.org/roadmap/, Vitalik Buterin, "Possible futures of the Ethereum protocol" series

---

## 2. Fork history and timeline

| Fork | Date | Status | Key contents |
|---|---|---|---|
| Paris (The Merge) | 2022-09-15 | Live | PoS transition, Beacon Chain merge, difficulty bomb removal |
| Shapella | 2023-04-12 | Live | Staking withdrawals (EIP-4895), warm COINBASE (EIP-3651) |
| Dencun | 2024-03-13 | Live | Proto-Danksharding (EIP-4844, blob transactions), transient storage (EIP-1153), beacon block root (EIP-4788) |
| Pectra | 2025-05-07 | Live | EOAs gain smart-contract functionality (EIP-7702), higher max effective balance (validator consolidation), blob target 3 to 6 (max 9) |
| **Fusaka** | 2025-12-03 | Live | **PeerDAS (EIP-7594)**, BPO forks (EIP-7892), gas-limit target raised toward 60M, transaction gas cap (EIP-7825), and more |
| **Glamsterdam** | Targeted Q4 2026 | In development (devnet) | **ePBS (EIP-7732)**, **BALs (EIP-7928)**, state gas repricing (EIP-8037/8038), churn-limit expansion (EIP-8061), and more |
| Hegota | Targeted 2027 | Planning stage | **FOCIL (EIP-7805)**, EIP-8141 (native account abstraction) under consideration |

Note: Glamsterdam's name combines the execution-layer upgrade name "Amsterdam" (a past Devconnect location) with the consensus-layer name "Gloas" (a star name). Hegota follows the same EL/CL name-combination convention.

**Sources**: ethereum.org/roadmap/, ethereum.org/roadmap/fusaka/, ethereum.org/roadmap/glamsterdam/

---

## 3. Fusaka (live, 2025-12-03) in detail

### Scale Blobs
- **PeerDAS (EIP-7594)** — Fusaka's headline feature. With data availability sampling, each full node only needs to hold 1/8 of blob data, a theoretical 8x scale-up. Any 50% of the data suffices to reconstruct the rest (failure probability driven down to roughly 1e-20 to 1e-24).
- **Blob-Parameter-Only (BPO) forks (EIP-7892)** — a mechanism to raise blob target/max counts via lightweight mini-forks without waiting for a major upgrade. The blob target went from 3 (Dencun) to 6 (Pectra) and can now be raised independently.
- **Blob base-fee floor (EIP-7918)** — pins a minimum blob fee tied to execution cost, preventing the blob fee market from collapsing to 1 wei and losing its function as a price signal.

### Scale L1
- **History expiry (EIP-7642)** — execution clients began supporting partial history expiry in July 2025, allowing pre-Merge history to be dropped from local storage and reducing disk usage.
- **MODEXP input bound (EIP-7823)** — caps input size at 8192 bits (1024 bytes).
- **Transaction gas limit cap (EIP-7825)** — caps any single transaction at 2^24 (16,777,216) gas.
- **MODEXP gas repricing (EIP-7883)** — raises minimum gas and other constants to match real computational cost.
- **RLP block size cap (EIP-7934)** — 10 MiB total, with 2 MiB reserved for consensus data.
- **Default gas limit target of 60M (EIP-7935)** — the first systematic push to raise the default gas limit since the Merge (30M to 36M to 45M to a 60M target).

**Sources**: ethereum.org/roadmap/fusaka/, EIP-7594/7892/7918/7642/7823/7825/7883/7934/7935

---

## 4. Glamsterdam (in development, targeted Q4 2026) in detail

Three goals: (1) parallelization, (2) expanded capacity, (3) curbing database bloat (sustainability).

### Scale L1 and parallel processing
- **Enshrined Proposer-Builder Separation (ePBS, EIP-7732)** — a headliner. Enshrines the proposer/builder handoff directly in the protocol, removing the need for off-protocol middleware such as MEV-Boost-style relays. Extends the block propagation window from about 2 seconds to about 9 seconds. Introduces a Payload Timeliness Committee (PTC) and dual-deadline logic.
- **Block-Level Access Lists (BALs, EIP-7928)** — the other headliner. Lists every account/state access a block's transactions will touch, along with post-execution values, up front. Enables parallel disk reads and executionless sync (copying final results instead of replaying every transaction).
- **eth/71 Block Access List Exchange (EIP-8159)** — the required networking companion that lets nodes actually exchange BALs over the peer-to-peer network.

### Network sustainability
- **State creation gas cost increase (EIP-8037)** — ties the cost of creating new accounts/contracts to their actual long-term storage burden. Introduces a dedicated state-gas "reservoir" model targeting a predictable state-growth rate of about 120 GiB/year.
- **State-access gas cost update (EIP-8038)** — corrects the underpricing of state-access opcodes such as `EXTCODESIZE`/`EXTCODECOPY`, improving DoS resistance.

### Network resilience
- **Exclude slashed validators from proposing (EIP-8045)**
- **Exit/consolidation churn expansion (EIP-8061)** — at current staking levels, roughly 4x more exit capacity and 2x more consolidation capacity; the weak subjectivity period shrinks from about 15.7 days to about 7 days.

### UX and developer experience
- **Intrinsic gas reduction (EIP-2780)** — cuts the base fee for a simple ETH transfer by up to 71%.
- **Deterministic Factory Predeploy (EIP-7997)** — places a deterministic factory contract at address 0x12 on every participating EVM chain, so smart-contract wallets can share the same address across L1 and L2s.

Additional EIPs under devnet testing: EIP-7610, 7688, 7778, 7843, 7976, 7981, 8024, 8246, 8282 (see Forkcast for the latest status).

**Sources**: ethereum.org/roadmap/glamsterdam/, EIP-7732/7928/8159/8037/8038/8045/8061/2780/7997

---

## 5. The future of the state tree: Verkle trees to a binary state tree (important shift)

Verkle trees (Banderwagon elliptic curve plus IPA commitments) were long considered the centerpiece of The Verge. That is now changing significantly:

- In March 2026, Vitalik published a proposal for a sweeping execution-layer overhaul centered on **EIP-7864 (binary state tree)**, which would replace the current hexary Keccak Merkle Patricia Trie with a binary tree built on a more efficient hash function (Blake3 or a Poseidon variant). Credit goes to Guillaume Ballet and others; the EIP has been in draft since January 2025.
- A binary tree produces Merkle branches about 4x shorter than Verkle trees, and switching the hash function to Blake3 or Poseidon could improve proving efficiency by a further 3x to 100x (though Poseidon still needs more security review).
- **The main driver of the shift is quantum resistance**: Verkle trees' Banderwagon/IPA construction is elliptic-curve based and therefore vulnerable to Shor's algorithm on a sufficiently large quantum computer. Interest in binary trees picked back up from mid-2024 onward largely because of this concern.
- Vitalik argues that the state tree and the VM together account for more than 80% of Ethereum's proving bottleneck, making both changes "basically mandatory."

This state-tree change also underpins stateless clients (verification via witnesses instead of a full local state database).

**Sources**: The Block (2026-03-01), EIP-7864, ethereum.org/roadmap/verkle-trees/, ethereum.org/roadmap/statelessness/

---

## 6. History and state expiry (The Purge)

- **History expiry (EIP-4444)** — allows nodes to drop pre-Merge block history from local storage. Partial support began in July 2025 and became mandatory in Fusaka. Dropped history is expected to be served by the Portal Network (a peer-to-peer network that distributes history across nodes, still in development) or by external archives (altruistic actors, DAOs, etc.). Community buy-in, more than the technology itself, is the main obstacle.
- **State expiry** — removes state that has not been touched in a long time from the active set (rent-based or time-based). Still in the research phase. If stateless clients and history expiry land first, the need for state expiry may diminish.
- **Address space extension** — a proposal to widen addresses from 20 to 32 bytes to store resurrection metadata needed for state expiry.

**Sources**: ethereum.org/roadmap/statelessness/, EIP-4444, EIP-7642

---

## 7. Account abstraction

| EIP / mechanism | Status | Description |
|---|---|---|
| **ERC-4337** | Live since 2023-03 | Account abstraction without protocol changes: UserOperation objects plus an EntryPoint contract. Over 26 million smart accounts and 170 million UserOperations processed to date. |
| **EIP-7702** | Live (Pectra, 2025-05-07) | Lets an EOA be represented by the code of an existing smart contract, gaining batching, gas sponsorship, and recovery mechanisms while remaining an EOA. |
| **EIP-8141** (native account abstraction) | Hegota candidate (targeted 2027) | Protocol-native account abstraction. **Key to the quantum-resistance migration path**: gives individual accounts "signature agility," letting each account switch to a post-quantum signature scheme without waiting for a protocol-wide migration. |

**Sources**: ethereum.org/roadmap/account-abstraction/, ethereum.org/roadmap/security/quantum-resistance/, EIP-7702/4337/8141

---

## 8. Censorship resistance and MEV (The Scourge)

- **FOCIL — Fork-Choice enforced Inclusion Lists (EIP-7805)** — the consensus-layer headliner for Hegota. A committee-based design that lets many validators force specific transactions into a block, guarding against builder-level censorship. Proposed in June 2024, building on the forward-inclusion-list work from the Pectra era.
- ePBS (EIP-7732, described above) also reduces reliance on MEV-related middleware such as MEV-Boost.

**Sources**: ethereum.org/roadmap/security/, EIP-7805, ethresear.ch FOCIL thread

---

## 9. Quantum resistance and Lean Ethereum (important area)

In February 2026, Vitalik published four areas of Ethereum's cryptography that are vulnerable to a sufficiently powerful quantum computer, along with a response plan. The Ethereum Foundation formed a dedicated Post-Quantum Security team in January 2026 (led by Thomas Coratger; progress tracked at pq.ethereum.org).

| Vulnerable area | Why it is vulnerable | Response plan |
|---|---|---|
| **Consensus-layer BLS signatures** | Elliptic-curve pairings are broken by Shor's algorithm | **leanXMSS** (a hash-based, quantum-safe signature scheme) plus **leanVM** (a SNARK-based aggregation engine that compresses signatures roughly 250x, offsetting XMSS's much larger size — about 3000 bytes versus BLS's 96 bytes) |
| **KZG commitments for data availability** | Rely on elliptic-curve pairings | Safe today thanks to the trusted-setup design, even against later quantum attacks. Long term: migrate to STARK-based or lattice-based commitments, both still under research |
| **EOA signatures (ECDSA)** | Any account that has sent a transaction has its public key exposed on-chain, creating "harvest now, decrypt later" risk | Rather than one protocol-wide cutover, **EIP-8141's signature agility** lets accounts migrate individually |
| **Application-layer ZK proofs** | Many SNARKs rely on elliptic-curve pairings | STARKs (hash-based, already quantum-resistant) are seeing natural ecosystem adoption |

Mapping to the finalized NIST standards (August 2024):

| Standard | Name | Type | Use case |
|---|---|---|---|
| FIPS 203 | ML-KEM | Lattice-based | Key encapsulation |
| FIPS 204 | ML-DSA (Dilithium) | Lattice-based | Digital signatures |
| FIPS 205 | SLH-DSA (SPHINCS+) | Hash-based | Digital signatures |

Provisional migration milestones (names and order may change): I* (PQ key registry) leads to J* (PQ signature-verification precompiles) leads to L* (PQ attestations and real-time consensus-layer proofs via leanVM) leads to M* (full PQ signature aggregation and PQ-safe blob commitments). Core infrastructure is targeted for completion around 2029.

leanXMSS, leanVM, leanSpec (Python), leanSig (Rust), and leanMultisig are all open source under the `leanEthereum` GitHub organization. More than ten client teams participate in weekly PQ interoperability devnets, including Lighthouse, Grandine, Zeam, Ream Labs, and PierTwo.

**Sources**: ethereum.org/roadmap/security/quantum-resistance/, pq.ethereum.org, NIST PQC standards

---

## 10. The future of the execution layer: RISC-V versus WASM

Since April 2025, Vitalik has floated replacing the EVM with the RISC-V instruction set, arguing that many ZK provers already use RISC-V internally. His three-stage rollout: (1) RISC-V only for precompiles, (2) users can deploy RISC-V contracts, (3) the EVM itself is eventually implemented as a RISC-V smart contract and treated as legacy.

In November 2025, Offchain Labs (the team behind Arbitrum) published a rebuttal arguing WASM is the better long-term choice, on the grounds that the "delivery ISA" and the "proving ISA" do not need to be the same thing. This debate has not reached broad consensus and remains considerably more speculative than the state-tree change (EIP-7864).

**Sources**: The Block (2026-03-01, 2025-04-20, 2025-11-23)

---

## 11. EIP quick reference

| EIP | Description | Target fork | Status |
|---|---|---|---|
| 4844 | Proto-Danksharding (blob transactions) | Dencun | Live |
| 7702 | EOAs gain smart-contract functionality | Pectra | Live |
| 7594 | PeerDAS | Fusaka | Live |
| 7892 | Blob-Parameter-Only forks | Fusaka | Live |
| 7642 | History expiry | Fusaka (partial support from 2025-07) | Live |
| 7825 | Transaction gas limit cap | Fusaka | Live |
| 7883 | MODEXP gas repricing | Fusaka | Live |
| 7934 | RLP block size cap | Fusaka | Live |
| 7935 | Default gas limit target of 60M | Fusaka | Live |
| 7732 | Enshrined PBS (ePBS) | Glamsterdam | Devnet |
| 7928 | Block-Level Access Lists (BALs) | Glamsterdam | Devnet |
| 8159 | eth/71 Block Access List Exchange | Glamsterdam | Devnet |
| 8037 | State creation gas cost increase | Glamsterdam | Devnet |
| 8038 | State-access gas cost update | Glamsterdam | Devnet |
| 8045 | Exclude slashed validators from proposing | Glamsterdam | Devnet |
| 8061 | Exit/consolidation churn expansion | Glamsterdam | Devnet |
| 2780 | Intrinsic gas reduction | Glamsterdam | Devnet |
| 7997 | Deterministic Factory Predeploy | Glamsterdam | Devnet |
| 7864 | Binary state tree (Verkle tree alternative) | Undecided (post-Glamsterdam or Hegota) | Draft (since 2025-01) |
| 7805 | FOCIL | Hegota | Planned (headliner) |
| 8141 | Native account abstraction (key to quantum migration) | Hegota candidate | Under consideration |
| 4444 | History expiry (full specification) | Undecided | Under discussion (community buy-in is the harder problem, not the tech) |
| — | State expiry | Undecided | Research phase |
| — | RISC-V VM replacement | Undecided (long term) | Speculative, no consensus |

---

## 12. Sources

- https://ethereum.org/roadmap/
- https://ethereum.org/roadmap/fusaka/
- https://ethereum.org/roadmap/glamsterdam/
- https://ethereum.org/roadmap/scaling/
- https://ethereum.org/roadmap/statelessness/
- https://ethereum.org/roadmap/verkle-trees/
- https://ethereum.org/roadmap/account-abstraction/
- https://ethereum.org/roadmap/security/
- https://ethereum.org/roadmap/security/quantum-resistance/
- https://eips.ethereum.org/EIPS/{7864, 7732, 7928, 8159, 8037, 8038, 8045, 8061, 2780, 7997, 7805, 8141, 4444, 7594, 7892, 7642, 7825, 7883, 7934, 7935, 7918, 7823}
- https://www.theblock.co/news/ecosystems/2026-03-01-vitalik-buterin-lays-out-a-two-part-plan-to-overhaul-ethereums-execution-layer-from-the-ground-up-391681
- https://x.com/VitalikButerin/status/2028158949720252574 (binary tree / RISC-V proposal)
- https://x.com/VitalikButerin/status/2027075026378543132 (four areas of quantum vulnerability)
- https://pq.ethereum.org/
- https://forkcast.org/upgrade/{fusaka,glamsterdam}
- https://ethresear.ch/t/fork-choice-enforced-inclusion-lists-focil-a-simple-committee-based-inclusion-list-proposal/19870
