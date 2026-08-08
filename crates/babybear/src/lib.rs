use field_core::Fp;

pub use field_core::FpBackend;

// The BabyBear prime, p = 2^31 - 2^27 + 1 = 15*2^27 + 1. Small enough that
// a*b for a, b < p fits in 62 bits, so plain u64 arithmetic (no montgomery
// tricks, no bigint) is enough for every operation here.
const BABYBEAR_P: u64 = 0x7800_0001;

pub struct BabyBearBackend;

impl FpBackend for BabyBearBackend {
    type Repr = u64;

    const MODULUS: u64 = BABYBEAR_P;

    fn add(a: u64, b: u64) -> u64 {
        let sum = a + b;

        if sum >= Self::MODULUS { sum - Self::MODULUS } else { sum }
    }

    fn sub(a: u64, b: u64) -> u64 {
        if a >= b { a - b } else { a + Self::MODULUS - b }
    }

    fn mul(a: u64, b: u64) -> u64 {
        (a * b) % Self::MODULUS
    }

    fn neg(a: u64) -> u64 {
        if a == 0 { 0 } else { Self::MODULUS - a }
    }

    // Fermat's little theorem: a^(p-2) == a^-1 (mod p).
    fn inverse(a: u64) -> u64 {
        assert!(a != 0, "cannot invert zero in a field");

        let mut result = 1u64;
        let mut base = a;
        let mut e = Self::MODULUS - 2;

        while e != 0 {
            if e & 1 == 1 {
                result = Self::mul(result, base);
            }
            base = Self::square(base);
            e >>= 1;
        }

        result
    }

    fn one() -> u64 {
        1
    }
}

pub type BabyBear = Fp<BabyBearBackend>;

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn babybear_modulus_is_correct() {
        assert_eq!(BabyBearBackend::MODULUS, 2_013_265_921);
    }

    fn fb(v: u64) -> BabyBear {
        BabyBear::new(v % BabyBearBackend::MODULUS)
    }

    #[test]
    fn babybear_mul_max_values() {
        let max = BabyBearBackend::MODULUS - 1;
        let expected = ((max as u128 * max as u128) % BabyBearBackend::MODULUS as u128) as u64;

        assert_eq!((fb(max) * fb(max)).value, expected);
    }

    proptest! {
        #[test]
        fn babybear_add_matches_u128_mod(a in any::<u64>(), b in any::<u64>()) {
            let expected = ((fb(a).value as u128 + fb(b).value as u128) % BabyBearBackend::MODULUS as u128) as u64;

            prop_assert_eq!((fb(a) + fb(b)).value, expected);
        }

        #[test]
        fn babybear_sub_matches_u128_mod(a in any::<u64>(), b in any::<u64>()) {
            let p = BabyBearBackend::MODULUS as u128;
            let expected = ((fb(a).value as u128 + p - fb(b).value as u128) % p) as u64;

            prop_assert_eq!((fb(a) - fb(b)).value, expected);
        }

        #[test]
        fn babybear_mul_matches_u128_mod(a in any::<u64>(), b in any::<u64>()) {
            let expected = ((fb(a).value as u128 * fb(b).value as u128) % BabyBearBackend::MODULUS as u128) as u64;

            prop_assert_eq!((fb(a) * fb(b)).value, expected);
        }

        #[test]
        fn babybear_neg_matches_u128_mod(a in any::<u64>()) {
            let p = BabyBearBackend::MODULUS as u128;
            let expected = ((p - fb(a).value as u128) % p) as u64;

            prop_assert_eq!((-fb(a)).value, expected);
        }

        #[test]
        fn babybear_square_matches_mul(a in any::<u64>()) {
            prop_assert_eq!(fb(a).square().value, (fb(a) * fb(a)).value);
        }

        #[test]
        fn babybear_result_is_canonical(a in any::<u64>(), b in any::<u64>()) {
            prop_assert!((fb(a) + fb(b)).value < BabyBearBackend::MODULUS);
            prop_assert!((fb(a) - fb(b)).value < BabyBearBackend::MODULUS);
            prop_assert!((fb(a) * fb(b)).value < BabyBearBackend::MODULUS);
        }

        #[test]
        fn babybear_inverse_times_self_is_one(a in 1..BabyBearBackend::MODULUS) {
            prop_assert_eq!((fb(a) * fb(a).inverse()).value, 1);
        }

        #[test]
        fn babybear_div_is_inverse_of_mul(a in any::<u64>(), b in 1..BabyBearBackend::MODULUS) {
            prop_assert_eq!(((fb(a) * fb(b)) / fb(b)).value, fb(a).value);
        }

        #[test]
        fn babybear_pow_matches_repeated_mul(a in any::<u64>(), exp in 0u32..12) {
            let mut expected = fb(1);
            for _ in 0..exp {
                expected = expected * fb(a);
            }

            prop_assert_eq!(fb(a).pow(exp as u64).value, expected.value);
        }

        #[test]
        fn babybear_fermats_little_theorem(a in 1..BabyBearBackend::MODULUS) {
            // a^(p-1) == 1 (mod p) for prime p and a not divisible by p.
            prop_assert_eq!(fb(a).pow(BabyBearBackend::MODULUS - 1).value, 1);
        }
    }
}
