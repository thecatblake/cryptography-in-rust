use field_core::{Fp, WideArithmeticBackend, WideFieldConfig};

pub use field_core::FpBackend;

// The Mersenne31 prime, p = 2^31 - 1. Being a Mersenne prime means
// 2^31 == 1 (mod p) by construction (p + 1 = 2^31), so a product's full
// 62-bit width folds down to canonical range with masks/shifts/adds only --
// no division, and no Montgomery-style multiply-by-N_PRIME either. Unlike
// Goldilocks (whose fold needs 128-bit intermediates) or BabyBear (whose
// weaker 2-adic relation only trims 4 bits per fold), M31's relation is a
// single clean 31/31 split that fully collapses in two folds.
const MERSENNE31_P: u32 = (1 << 31) - 1;

// x < 2^62 (the product of two values < p < 2^31) splits into
// x_hi*2^31 + x_lo, which is congruent to x_hi + x_lo (mod p). That sum is
// < 2^32, so a second identical fold (this time against a 32-bit value)
// finishes the reduction, leaving a result in [0, p+1] that one conditional
// subtract makes canonical.
fn mersenne31_reduce62(x: u64) -> u32 {
    let x_lo = (x & MERSENNE31_P as u64) as u32;
    let x_hi = (x >> 31) as u32;
    let sum = x_lo + x_hi;

    let sum_lo = sum & MERSENNE31_P;
    let sum_hi = sum >> 31;
    let result = sum_lo + sum_hi;

    if result >= MERSENNE31_P { result - MERSENNE31_P } else { result }
}

pub struct Mersenne31Config;

impl WideFieldConfig for Mersenne31Config {
    type Repr = u32;

    const MODULUS: u32 = MERSENNE31_P;

    fn mul(a: u32, b: u32) -> u32 {
        mersenne31_reduce62(a as u64 * b as u64)
    }
}

pub type Mersenne31Backend = WideArithmeticBackend<Mersenne31Config>;
pub type Mersenne31 = Fp<Mersenne31Backend>;

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn mersenne31_modulus_is_correct() {
        assert_eq!(Mersenne31Backend::MODULUS, 2_147_483_647);
    }

    fn fm(v: u32) -> Mersenne31 {
        Mersenne31::new(v % Mersenne31Backend::MODULUS)
    }

    #[test]
    fn mersenne31_mul_max_values() {
        let max = Mersenne31Backend::MODULUS - 1;
        let expected = ((max as u64 * max as u64) % Mersenne31Backend::MODULUS as u64) as u32;

        assert_eq!((fm(max) * fm(max)).value, expected);
    }

    #[test]
    fn mersenne31_reduce62_matches_u64_mod_at_boundaries() {
        let p = MERSENNE31_P as u64;
        // 0, p*p (should fold to exactly 0), and the max possible 62-bit
        // product (p-1)*(p-1), which forces both folds and the final
        // conditional subtract.
        for x in [0u64, p * p, (p - 1) * (p - 1), 1u64 << 61, (1u64 << 62) - 1] {
            assert_eq!(mersenne31_reduce62(x) as u64, x % p);
        }
    }

    proptest! {
        #[test]
        fn mersenne31_reduce62_matches_u64_mod(a in 0..MERSENNE31_P, b in 0..MERSENNE31_P) {
            let x = a as u64 * b as u64;

            prop_assert_eq!(mersenne31_reduce62(x) as u64, x % MERSENNE31_P as u64);
        }

        #[test]
        fn mersenne31_add_matches_u64_mod(a in any::<u32>(), b in any::<u32>()) {
            let expected = ((fm(a).value as u64 + fm(b).value as u64) % Mersenne31Backend::MODULUS as u64) as u32;

            prop_assert_eq!((fm(a) + fm(b)).value, expected);
        }

        #[test]
        fn mersenne31_sub_matches_u64_mod(a in any::<u32>(), b in any::<u32>()) {
            let p = Mersenne31Backend::MODULUS as u64;
            let expected = ((fm(a).value as u64 + p - fm(b).value as u64) % p) as u32;

            prop_assert_eq!((fm(a) - fm(b)).value, expected);
        }

        #[test]
        fn mersenne31_mul_matches_u64_mod(a in any::<u32>(), b in any::<u32>()) {
            let expected = ((fm(a).value as u64 * fm(b).value as u64) % Mersenne31Backend::MODULUS as u64) as u32;

            prop_assert_eq!((fm(a) * fm(b)).value, expected);
        }

        #[test]
        fn mersenne31_neg_matches_u64_mod(a in any::<u32>()) {
            let p = Mersenne31Backend::MODULUS as u64;
            let expected = ((p - fm(a).value as u64) % p) as u32;

            prop_assert_eq!((-fm(a)).value, expected);
        }

        #[test]
        fn mersenne31_square_matches_mul(a in any::<u32>()) {
            prop_assert_eq!(fm(a).square().value, (fm(a) * fm(a)).value);
        }

        #[test]
        fn mersenne31_result_is_canonical(a in any::<u32>(), b in any::<u32>()) {
            prop_assert!((fm(a) + fm(b)).value < Mersenne31Backend::MODULUS);
            prop_assert!((fm(a) - fm(b)).value < Mersenne31Backend::MODULUS);
            prop_assert!((fm(a) * fm(b)).value < Mersenne31Backend::MODULUS);
        }

        #[test]
        fn mersenne31_inverse_times_self_is_one(a in 1..Mersenne31Backend::MODULUS) {
            prop_assert_eq!((fm(a) * fm(a).inverse()).value, 1);
        }

        #[test]
        fn mersenne31_div_is_inverse_of_mul(a in any::<u32>(), b in 1..Mersenne31Backend::MODULUS) {
            prop_assert_eq!(((fm(a) * fm(b)) / fm(b)).value, fm(a).value);
        }

        #[test]
        fn mersenne31_pow_matches_repeated_mul(a in any::<u32>(), exp in 0u32..12) {
            let mut expected = fm(1);
            for _ in 0..exp {
                expected = expected * fm(a);
            }

            prop_assert_eq!(fm(a).pow(exp).value, expected.value);
        }

        #[test]
        fn mersenne31_fermats_little_theorem(a in 1..Mersenne31Backend::MODULUS) {
            // a^(p-1) == 1 (mod p) for prime p and a not divisible by p.
            prop_assert_eq!(fm(a).pow(Mersenne31Backend::MODULUS - 1).value, 1);
        }
    }
}
