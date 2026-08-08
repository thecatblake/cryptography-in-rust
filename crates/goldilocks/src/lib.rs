use field_core::{Fp, NativeArithmeticBackend, NativeFieldConfig};

pub use field_core::FpBackend;

// The Goldilocks prime, p = 2^64 - 2^32 + 1. Fits in a single u64 limb, so
// this backend skips bigint entirely and does arithmetic natively.
const GOLDILOCKS_P: u64 = 0xFFFF_FFFF_0000_0001;
// 2^64 mod p: p = 2^64 - 2^32 + 1, so 2^64 == 2^32 - 1 (mod p).
const GOLDILOCKS_EPSILON: u64 = 0xFFFF_FFFF;

// x = x_hi*2^64 + x_lo == x_hi*EPSILON + x_lo (mod p). Splitting
// x_hi = x_hi_hi*2^32 + x_hi_lo and using 2^32*EPSILON == -1 (mod p) reduces
// that further to x_lo - x_hi_hi + x_hi_lo*EPSILON, where every intermediate
// term fits in a u64 (x_hi_lo, EPSILON < 2^32, so their product is < 2^64).
fn goldilocks_reduce128(x: u128) -> u64 {
    let x_lo = x as u64;
    let x_hi = (x >> 64) as u64;

    let x_hi_hi = x_hi >> 32;
    let x_hi_lo = x_hi & 0xFFFF_FFFF;

    let (t0, borrow) = x_lo.overflowing_sub(x_hi_hi);
    let t0 = if borrow { t0.wrapping_sub(GOLDILOCKS_EPSILON) } else { t0 };

    let t1 = x_hi_lo * GOLDILOCKS_EPSILON;

    let (t2, carry) = t0.overflowing_add(t1);
    let t2 = if carry { t2.wrapping_add(GOLDILOCKS_EPSILON) } else { t2 };

    if t2 >= GOLDILOCKS_P { t2 - GOLDILOCKS_P } else { t2 }
}

pub struct GoldilocksConfig;

impl NativeFieldConfig for GoldilocksConfig {
    type Repr = u64;

    const MODULUS: u64 = GOLDILOCKS_P;

    fn mul(a: u64, b: u64) -> u64 {
        goldilocks_reduce128(a as u128 * b as u128)
    }
}

pub type GoldilocksBackend = NativeArithmeticBackend<GoldilocksConfig>;
pub type Goldilocks = Fp<GoldilocksBackend>;

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn goldilocks_modulus_is_correct() {
        assert_eq!(GoldilocksBackend::MODULUS, 0xFFFF_FFFF_0000_0001);
    }

    fn fg(v: u64) -> Goldilocks {
        Goldilocks::new(v % GoldilocksBackend::MODULUS)
    }

    #[test]
    fn goldilocks_add_overflows_native_u64() {
        // (p-1) + (p-1) overflows a u64 add outright; exercises the carry path.
        let expected = GoldilocksBackend::MODULUS - 2;

        assert_eq!((fg(GoldilocksBackend::MODULUS - 1) + fg(GoldilocksBackend::MODULUS - 1)).value, expected);
    }

    #[test]
    fn goldilocks_mul_max_values() {
        let max = GoldilocksBackend::MODULUS - 1;
        let expected = ((max as u128 * max as u128) % GoldilocksBackend::MODULUS as u128) as u64;

        assert_eq!((fg(max) * fg(max)).value, expected);
    }

    proptest! {
        #[test]
        fn goldilocks_add_matches_u128_mod(a in any::<u64>(), b in any::<u64>()) {
            let expected = ((fg(a).value as u128 + fg(b).value as u128) % GoldilocksBackend::MODULUS as u128) as u64;

            prop_assert_eq!((fg(a) + fg(b)).value, expected);
        }

        #[test]
        fn goldilocks_sub_matches_u128_mod(a in any::<u64>(), b in any::<u64>()) {
            let p = GoldilocksBackend::MODULUS as u128;
            let expected = ((fg(a).value as u128 + p - fg(b).value as u128) % p) as u64;

            prop_assert_eq!((fg(a) - fg(b)).value, expected);
        }

        #[test]
        fn goldilocks_mul_matches_u128_mod(a in any::<u64>(), b in any::<u64>()) {
            let expected = ((fg(a).value as u128 * fg(b).value as u128) % GoldilocksBackend::MODULUS as u128) as u64;

            prop_assert_eq!((fg(a) * fg(b)).value, expected);
        }

        #[test]
        fn goldilocks_neg_matches_u128_mod(a in any::<u64>()) {
            let p = GoldilocksBackend::MODULUS as u128;
            let expected = ((p - fg(a).value as u128) % p) as u64;

            prop_assert_eq!((-fg(a)).value, expected);
        }

        #[test]
        fn goldilocks_square_matches_mul(a in any::<u64>()) {
            prop_assert_eq!(fg(a).square().value, (fg(a) * fg(a)).value);
        }

        #[test]
        fn goldilocks_result_is_canonical(a in any::<u64>(), b in any::<u64>()) {
            prop_assert!((fg(a) + fg(b)).value < GoldilocksBackend::MODULUS);
            prop_assert!((fg(a) - fg(b)).value < GoldilocksBackend::MODULUS);
            prop_assert!((fg(a) * fg(b)).value < GoldilocksBackend::MODULUS);
        }

        #[test]
        fn goldilocks_inverse_times_self_is_one(a in 1..GoldilocksBackend::MODULUS) {
            prop_assert_eq!((fg(a) * fg(a).inverse()).value, 1);
        }

        #[test]
        fn goldilocks_div_is_inverse_of_mul(a in any::<u64>(), b in 1..GoldilocksBackend::MODULUS) {
            prop_assert_eq!(((fg(a) * fg(b)) / fg(b)).value, fg(a).value);
        }

        #[test]
        fn goldilocks_pow_matches_repeated_mul(a in any::<u64>(), exp in 0u32..12) {
            let mut expected = fg(1);
            for _ in 0..exp {
                expected = expected * fg(a);
            }

            prop_assert_eq!(fg(a).pow(exp as u64).value, expected.value);
        }

        #[test]
        fn goldilocks_fermats_little_theorem(a in 1..GoldilocksBackend::MODULUS) {
            // a^(p-1) == 1 (mod p) for prime p and a not divisible by p.
            prop_assert_eq!(fg(a).pow(GoldilocksBackend::MODULUS - 1).value, 1);
        }
    }
}
