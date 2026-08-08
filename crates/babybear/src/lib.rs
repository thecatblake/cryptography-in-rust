use field_core::{Fp, NativeArithmeticBackend, NativeFieldConfig};

pub use field_core::FpBackend;

// The BabyBear prime, p = 2^31 - 2^27 + 1 = 15*2^27 + 1. Fits in 31 bits, so
// (unlike Goldilocks, whose prime is nearly u64::MAX) it's represented as a
// u32 here to match how real BabyBear implementations pack it (e.g. 4 lanes
// per 128-bit SIMD register instead of 2). a*b for a, b < p fits in 62 bits,
// so widening to u64 for the multiplication reduction is enough.
const BABYBEAR_P: u32 = 0x7800_0001;

pub struct BabyBearConfig;

impl NativeFieldConfig for BabyBearConfig {
    type Repr = u32;

    const MODULUS: u32 = BABYBEAR_P;

    fn mul(a: u32, b: u32) -> u32 {
        ((a as u64 * b as u64) % Self::MODULUS as u64) as u32
    }
}

pub type BabyBearBackend = NativeArithmeticBackend<BabyBearConfig>;
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
        BabyBear::new(v % BabyBearBackend::MODULUS)
    }

    #[test]
    fn babybear_mul_max_values() {
        let max = BabyBearBackend::MODULUS - 1;
        let expected = ((max as u64 * max as u64) % BabyBearBackend::MODULUS as u64) as u32;

        assert_eq!((fb(max) * fb(max)).value, expected);
    }

    proptest! {
        #[test]
        fn babybear_add_matches_u64_mod(a in any::<u32>(), b in any::<u32>()) {
            let expected = ((fb(a).value as u64 + fb(b).value as u64) % BabyBearBackend::MODULUS as u64) as u32;

            prop_assert_eq!((fb(a) + fb(b)).value, expected);
        }

        #[test]
        fn babybear_sub_matches_u64_mod(a in any::<u32>(), b in any::<u32>()) {
            let p = BabyBearBackend::MODULUS as u64;
            let expected = ((fb(a).value as u64 + p - fb(b).value as u64) % p) as u32;

            prop_assert_eq!((fb(a) - fb(b)).value, expected);
        }

        #[test]
        fn babybear_mul_matches_u64_mod(a in any::<u32>(), b in any::<u32>()) {
            let expected = ((fb(a).value as u64 * fb(b).value as u64) % BabyBearBackend::MODULUS as u64) as u32;

            prop_assert_eq!((fb(a) * fb(b)).value, expected);
        }

        #[test]
        fn babybear_neg_matches_u64_mod(a in any::<u32>()) {
            let p = BabyBearBackend::MODULUS as u64;
            let expected = ((p - fb(a).value as u64) % p) as u32;

            prop_assert_eq!((-fb(a)).value, expected);
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
            prop_assert_eq!((fb(a) * fb(a).inverse()).value, 1);
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
            prop_assert_eq!(fb(a).pow(BabyBearBackend::MODULUS - 1).value, 1);
        }
    }
}
