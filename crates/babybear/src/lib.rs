use field_core::{Fp, MontFieldConfig, MontWideBackend};

pub use field_core::{Field, FpBackend, Frobenius};

mod extension;
pub use extension::{BabyBearFp2, BabyBearFp6, BabyBearFp12};

// The BabyBear prime, p = 2^31 - 2^27 + 1 = 15*2^27 + 1. Represented as u32
// (unlike Goldilocks' u64, since p fits in 31 bits) to match how real
// BabyBear implementations pack it -- e.g. 4 lanes per 128-bit SIMD
// register instead of 2.
//
// Unlike Goldilocks, p doesn't get a cheap Solinas-style bit-trick
// reduction: Goldilocks' defining relation (2^64 == 2^32 - 1 mod p) halves
// a product's bit-width in a single fold, because its prime is shaped
// specifically for that. BabyBear's relation (2^31 == 2^27 - 1 mod p) only
// shrinks a product by 4 bits per fold (p was chosen for its 2-adicity,
// i.e. 2^27 | p-1, for NTTs -- not for reduction-friendliness), which would
// need ~8 folds to be worth anything. Montgomery form (R = 2^32) reduces in
// one pass regardless of the prime's shape, at the cost of representing
// values pre-multiplied by R.
const BABYBEAR_P: u32 = 0x7800_0001;

pub struct BabyBearConfig;

impl MontFieldConfig for BabyBearConfig {
    type Repr = u32;

    const MODULUS: u32 = BABYBEAR_P;
    // R2 = R^2 mod p and N_PRIME = -p^-1 mod R (R = 2^32), computed and
    // round-trip verified offline; see babybear_mont_roundtrip_is_identity.
    const R2: u32 = 0x45dd_dde3;
    const N_PRIME: u32 = 0x77ff_ffff;
}

pub type BabyBearBackend = MontWideBackend<BabyBearConfig>;
pub type BabyBear = Fp<BabyBearBackend>;

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn babybear_modulus_is_correct() {
        assert_eq!(BabyBearBackend::MODULUS, 2_013_265_921);
    }

    fn fb(v: u32) -> BabyBear {
        BabyBear::new(BabyBearBackend::to_mont(v % BabyBearBackend::MODULUS))
    }

    fn canonical(x: BabyBear) -> u32 {
        BabyBearBackend::from_mont(x.value)
    }

    #[test]
    fn babybear_mont_roundtrip_is_identity() {
        let v = 123_456_789 % BabyBearBackend::MODULUS;

        assert_eq!(BabyBearBackend::from_mont(BabyBearBackend::to_mont(v)), v);
    }

    #[test]
    fn babybear_mul_max_values() {
        let max = BabyBearBackend::MODULUS - 1;
        let expected = ((max as u64 * max as u64) % BabyBearBackend::MODULUS as u64) as u32;

        assert_eq!(canonical(fb(max) * fb(max)), expected);
    }

    proptest! {
        #[test]
        fn babybear_add_matches_u64_mod(a in any::<u32>(), b in any::<u32>()) {
            let a_r = a % BabyBearBackend::MODULUS;
            let b_r = b % BabyBearBackend::MODULUS;
            let expected = ((a_r as u64 + b_r as u64) % BabyBearBackend::MODULUS as u64) as u32;

            prop_assert_eq!(canonical(fb(a) + fb(b)), expected);
        }

        #[test]
        fn babybear_sub_matches_u64_mod(a in any::<u32>(), b in any::<u32>()) {
            let a_r = a % BabyBearBackend::MODULUS;
            let b_r = b % BabyBearBackend::MODULUS;
            let p = BabyBearBackend::MODULUS as u64;
            let expected = ((a_r as u64 + p - b_r as u64) % p) as u32;

            prop_assert_eq!(canonical(fb(a) - fb(b)), expected);
        }

        #[test]
        fn babybear_mul_matches_u64_mod(a in any::<u32>(), b in any::<u32>()) {
            let a_r = a % BabyBearBackend::MODULUS;
            let b_r = b % BabyBearBackend::MODULUS;
            let expected = ((a_r as u64 * b_r as u64) % BabyBearBackend::MODULUS as u64) as u32;

            prop_assert_eq!(canonical(fb(a) * fb(b)), expected);
        }

        #[test]
        fn babybear_neg_matches_u64_mod(a in any::<u32>()) {
            let a_r = a % BabyBearBackend::MODULUS;
            let p = BabyBearBackend::MODULUS as u64;
            let expected = ((p - a_r as u64) % p) as u32;

            prop_assert_eq!(canonical(-fb(a)), expected);
        }

        #[test]
        fn babybear_square_matches_mul(a in any::<u32>()) {
            prop_assert_eq!(fb(a).square().value, (fb(a) * fb(a)).value);
        }

        #[test]
        fn babybear_result_is_canonical(a in any::<u32>(), b in any::<u32>()) {
            prop_assert!((fb(a) + fb(b)).value < BabyBearBackend::MODULUS);
            prop_assert!((fb(a) - fb(b)).value < BabyBearBackend::MODULUS);
            prop_assert!((fb(a) * fb(b)).value < BabyBearBackend::MODULUS);
        }

        #[test]
        fn babybear_inverse_times_self_is_one(a in 1..BabyBearBackend::MODULUS) {
            prop_assert_eq!(canonical(fb(a) * fb(a).inverse()), 1);
        }

        #[test]
        fn babybear_div_is_inverse_of_mul(a in any::<u32>(), b in 1..BabyBearBackend::MODULUS) {
            prop_assert_eq!(((fb(a) * fb(b)) / fb(b)).value, fb(a).value);
        }

        #[test]
        fn babybear_pow_matches_repeated_mul(a in any::<u32>(), exp in 0u32..12) {
            let mut expected = fb(1);
            for _ in 0..exp {
                expected = expected * fb(a);
            }

            prop_assert_eq!(fb(a).pow(exp).value, expected.value);
        }

        #[test]
        fn babybear_fermats_little_theorem(a in 1..BabyBearBackend::MODULUS) {
            // a^(p-1) == 1 (mod p) for prime p and a not divisible by p.
            prop_assert_eq!(canonical(fb(a).pow(BabyBearBackend::MODULUS - 1)), 1);
        }
    }
}
