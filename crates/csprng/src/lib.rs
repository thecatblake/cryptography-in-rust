// rand_core-style byte/word source. there's no OS-entropy-unavailable case to surface here,
// unlike rand_core's split RngCore/TryRngCore.
pub trait Rand {
    fn next_u32(&mut self) -> u32;
    fn next_u64(&mut self) -> u64;
    fn next_bytes(&mut self, dest: &mut [u8]);
}

// Marker trait: opting in asserts the generator's output is
// indistinguishable from random to a computationally bounded adversary
// CryptoRng -- so callers that need randomness for key material can bound
// a generic parameter by `R: Rand + CryptoRng` and reject anything that
// hasn't made this claim.
pub trait CryptoRng: Rand {}

// Constructs a generator from explicit seed material. Two entry points,
// same split rand_core makes: from_seed for reproducible test vectors and
// protocol-level derivation (exact bytes matter), seed_from_u64 for quick
// throwaway seeding (only a u64 of entropy on hand, e.g. a test index).
pub trait SeedableRng: Sized {
    type Seed: Default + AsMut<[u8]>;

    fn from_seed(seed: Self::Seed) -> Self;
}
