use bigint::{U256, U512, Uint};
use bigint::math::mod_inv;
use field_core::FpBackend;
use std::marker::PhantomData;

pub use field_core::Fp;

// FpConfig associates a modulus with the type rather than each value.
// This makes field elements with different moduli different Rust types,
// and we can impose the operation under the same modulus at the compile time

// We need to separate the computation at compile time and runtime.
// Calculation methods are determined at compile time.
// The actual value is determined at runtime.
// So Fp holds a value and template type T for Fp holds calculation methods and template type for T holds constants.


pub trait FpConfig {
    const MODULUS: U256;
}

pub trait FpMontConfig {
    const MODULUS: U256;
    const R2: U256;
    const N_PRIME: U256;
}

// Scratch width for REDC's T + m*MODULUS step: both terms are < R*MODULUS < 2^512,
// so their sum can briefly need one bit more than U512 holds.
type Wide = Uint<9>;

pub struct DefaultBackend<T: FpConfig>(PhantomData<T>);
pub struct MontBackend<T: FpMontConfig>(PhantomData<T>);

impl<T: FpConfig> FpBackend for DefaultBackend<T> {
    type Repr = U256;

    const MODULUS: U256 = T::MODULUS;

    fn add(a: U256, b: U256) -> U256 {
        let mut sum = a + b;

        if sum >= Self::MODULUS {
            sum -= Self::MODULUS;
        }

        sum
    }

    fn sub(a: U256, b: U256) -> U256 {
        if a >= b {
            a - b
        } else {
            (a + Self::MODULUS) - b
        }
    }

    fn mul(a: U256, b: U256) -> U256 {
        let product: U512 = a * b;
        let modulus: U512 = Self::MODULUS.resize();

        (product % modulus).resize()
    }

    fn neg(a: U256) -> U256 {
        if a == U256::ZERO {
            U256::ZERO
        } else {
            Self::MODULUS - a
        }
    }

    fn inverse(a: U256) -> U256 {
        assert!(a != U256::ZERO, "cannot invert zero in a field");

        mod_inv(a, Self::MODULUS)
    }

    fn one() -> U256 {
        U256::ONE
    }

    fn square(a: U256) -> U256 {
        let product: U512 = a.square();
        let modulus: U512 = Self::MODULUS.resize();

        (product % modulus).resize()
    }
}

impl<T: FpMontConfig> MontBackend<T> {
    // REDC(t) = t * R^-1 mod MODULUS, computed without dividing by MODULUS.
    // Requires t < R * MODULUS, which holds for any t formed by multiplying
    // two values already reduced mod MODULUS.
    fn redc(t: U512) -> U256 {
        let t_low: U256 = t.resize();
        let m: U256 = (t_low * T::N_PRIME).resize();

        // t + m * MODULUS is guaranteed divisible by R (2^256) by construction of m.
        let sum: Wide = t.resize() + (m * T::MODULUS).resize();

        // Divide by R = 2^256 (i.e. drop the low 256 bits) via chained shifts,
        // since a single shift is limited to 63 bits at a time.
        let quotient: Wide = sum >> 63 >> 63 >> 63 >> 63 >> 4;

        let modulus_wide: Wide = T::MODULUS.resize();
        let reduced = if quotient >= modulus_wide {
            quotient - modulus_wide
        } else {
            quotient
        };

        reduced.resize()
    }

    fn to_mont(a: U256) -> U256 {
        Self::redc(a * T::R2)
    }

    fn from_mont(a: U256) -> U256 {
        Self::redc(a.resize())
    }
}

impl<T: FpMontConfig> FpBackend for MontBackend<T> {
    type Repr = U256;

    const MODULUS: U256 = T::MODULUS;

    fn add(a: U256, b: U256) -> U256 {
        let mut sum = a + b;

        if sum >= Self::MODULUS {
            sum -= Self::MODULUS;
        }

        sum
    }

    fn sub(a: U256, b: U256) -> U256 {
        if a >= b {
            a - b
        } else {
            (a + Self::MODULUS) - b
        }
    }

    fn mul(a: U256, b: U256) -> U256 {
        Self::redc(a * b)
    }

    fn neg(a: U256) -> U256 {
        if a == U256::ZERO {
            U256::ZERO
        } else {
            Self::MODULUS - a
        }
    }

    fn inverse(a: U256) -> U256 {
        assert!(a != U256::ZERO, "cannot invert zero in a field");

        let plain = Self::from_mont(a);
        let inv = mod_inv(plain, Self::MODULUS);

        Self::to_mont(inv)
    }

    fn one() -> U256 {
        Self::to_mont(U256::ONE)
    }

    fn square(a: U256) -> U256 {
        Self::redc(a.square())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    struct Mod17;
    impl FpConfig for Mod17 {
        const MODULUS: U256 = U256::from_u64(17);
    }

    type F17 = Fp<DefaultBackend<Mod17>>;

    fn fe(v: u64) -> F17 {
        F17::new(U256::from(v % 17))
    }

    #[test]
    fn add_no_reduction_needed() {
        let sum = fe(3) + fe(4);
        assert_eq!(sum.value, U256::from(7));
    }

    #[test]
    fn add_wraps_at_modulus() {
        let sum = fe(16) + fe(1);
        assert_eq!(sum.value, U256::from(0));
    }

    #[test]
    fn add_wraps_past_modulus() {
        let sum = fe(16) + fe(16);
        assert_eq!(sum.value, U256::from(15));
    }

    proptest! {
        #[test]
        fn add_result_is_canonical(a in 0u64..17, b in 0u64..17) {
            let sum = fe(a) + fe(b);
            prop_assert!(sum.value < Mod17::MODULUS);
        }

        #[test]
        fn add_matches_u64_mod(a in 0u64..17, b in 0u64..17) {
            let sum = fe(a) + fe(b);
            prop_assert_eq!(sum.value, U256::from((a + b) % 17));
        }

        #[test]
        fn add_zero_is_identity(a in 0u64..17) {
            let sum = fe(a) + fe(0);
            prop_assert_eq!(sum.value, U256::from(a));
        }

        #[test]
        fn add_commutative(a in 0u64..17, b in 0u64..17) {
            prop_assert_eq!((fe(a) + fe(b)).value, (fe(b) + fe(a)).value);
        }
    }

    #[test]
    fn sub_no_borrow_needed() {
        let diff = fe(7) - fe(3);
        assert_eq!(diff.value, U256::from(4));
    }

    #[test]
    fn sub_wraps_below_zero() {
        let diff = fe(3) - fe(7);
        assert_eq!(diff.value, U256::from(13));
    }

    proptest! {
        #[test]
        fn sub_result_is_canonical(a in 0u64..17, b in 0u64..17) {
            let diff = fe(a) - fe(b);
            prop_assert!(diff.value < Mod17::MODULUS);
        }

        #[test]
        fn sub_matches_u64_mod(a in 0u64..17, b in 0u64..17) {
            let diff = fe(a) - fe(b);
            let expected = (a + 17 - b) % 17;
            prop_assert_eq!(diff.value, U256::from(expected));
        }

        #[test]
        fn sub_is_inverse_of_add(a in 0u64..17, b in 0u64..17) {
            prop_assert_eq!(((fe(a) + fe(b)) - fe(b)).value, U256::from(a));
        }

        #[test]
        fn mul_result_is_canonical(a in 0u64..17, b in 0u64..17) {
            let product = fe(a) * fe(b);
            prop_assert!(product.value < Mod17::MODULUS);
        }

        #[test]
        fn mul_matches_u64_mod(a in 0u64..17, b in 0u64..17) {
            let product = fe(a) * fe(b);
            prop_assert_eq!(product.value, U256::from((a * b) % 17));
        }

        #[test]
        fn mul_zero_is_absorbing(a in 0u64..17) {
            prop_assert_eq!((fe(a) * fe(0)).value, U256::from(0));
        }

        #[test]
        fn mul_one_is_identity(a in 0u64..17) {
            prop_assert_eq!((fe(a) * fe(1)).value, U256::from(a));
        }

        #[test]
        fn mul_commutative(a in 0u64..17, b in 0u64..17) {
            prop_assert_eq!((fe(a) * fe(b)).value, (fe(b) * fe(a)).value);
        }

        #[test]
        fn square_matches_mul(a in 0u64..17) {
            prop_assert_eq!(fe(a).square().value, (fe(a) * fe(a)).value);
        }

        #[test]
        fn square_result_is_canonical(a in 0u64..17) {
            prop_assert!(fe(a).square().value < Mod17::MODULUS);
        }

        #[test]
        fn inverse_times_self_is_one(a in 1u64..17) {
            prop_assert_eq!((fe(a) * fe(a).inverse()).value, U256::from(1));
        }

        #[test]
        fn div_is_inverse_of_mul(a in 0u64..17, b in 1u64..17) {
            prop_assert_eq!(((fe(a) * fe(b)) / fe(b)).value, U256::from(a));
        }

        #[test]
        fn div_by_self_is_one(a in 1u64..17) {
            prop_assert_eq!((fe(a) / fe(a)).value, U256::from(1));
        }
    }

    #[test]
    #[should_panic(expected = "cannot invert zero")]
    fn inverse_of_zero_panics() {
        fe(0).inverse();
    }

    #[test]
    fn neg_zero_is_zero() {
        assert_eq!((-fe(0)).value, U256::from(0));
    }

    #[test]
    fn neg_nonzero_example() {
        assert_eq!((-fe(5)).value, U256::from(12));
    }

    proptest! {
        #[test]
        fn neg_result_is_canonical(a in 0u64..17) {
            prop_assert!((-fe(a)).value < Mod17::MODULUS);
        }

        #[test]
        fn neg_matches_u64_mod(a in 0u64..17) {
            let expected = (17 - a) % 17;
            prop_assert_eq!((-fe(a)).value, U256::from(expected));
        }

        #[test]
        fn add_neg_is_zero(a in 0u64..17) {
            prop_assert_eq!((fe(a) + (-fe(a))).value, U256::from(0));
        }

        #[test]
        fn double_neg_is_identity(a in 0u64..17) {
            prop_assert_eq!((-(-fe(a))).value, U256::from(a));
        }

        #[test]
        fn neg_matches_sub_from_zero(a in 0u64..17) {
            prop_assert_eq!((-fe(a)).value, (fe(0) - fe(a)).value);
        }
    }

    // If all the tests for DefaultBackend passes and the computation using the MontBackend and the one using DefaultBackend matches MontBackend passes all the tests intended.
    struct Mod17Mont;
    impl FpMontConfig for Mod17Mont {
        const MODULUS: U256 = U256::from_u64(17);
        // R = 2^256; ord(2 mod 17) = 8 and 8 | 256, so 2^256 ≡ 1 (mod 17) and R2 = 1.
        const R2: U256 = U256::from_u64(1);
        // -17^-1 mod 2^256, via the 2-adic Newton iteration (mirrors the classic
        // "17 * 0x0F...0F ≡ -1" byte-wise pattern extended to all 256 bits).
        const N_PRIME: U256 = U256::from_limbs([
            0x0F0F0F0F0F0F0F0F,
            0x0F0F0F0F0F0F0F0F,
            0x0F0F0F0F0F0F0F0F,
            0x0F0F0F0F0F0F0F0F,
        ]);
    }

    type MB17 = MontBackend<Mod17Mont>;
    type F17Mont = Fp<MB17>;

    fn fe_mont(v: u64) -> F17Mont {
        F17Mont::new(MB17::to_mont(U256::from(v % 17)))
    }

    fn canonical(x: F17Mont) -> U256 {
        MB17::from_mont(x.value)
    }

    proptest! {
        #[test]
        fn mont_roundtrip_is_identity(a in 0u64..17) {
            let x = U256::from(a);
            prop_assert_eq!(MB17::from_mont(MB17::to_mont(x)), x);
        }

        #[test]
        fn mont_add_matches_default(a in 0u64..17, b in 0u64..17) {
            prop_assert_eq!(canonical(fe_mont(a) + fe_mont(b)), (fe(a) + fe(b)).value);
        }

        #[test]
        fn mont_sub_matches_default(a in 0u64..17, b in 0u64..17) {
            prop_assert_eq!(canonical(fe_mont(a) - fe_mont(b)), (fe(a) - fe(b)).value);
        }

        #[test]
        fn mont_mul_matches_default(a in 0u64..17, b in 0u64..17) {
            prop_assert_eq!(canonical(fe_mont(a) * fe_mont(b)), (fe(a) * fe(b)).value);
        }

        #[test]
        fn mont_square_matches_default(a in 0u64..17) {
            prop_assert_eq!(canonical(fe_mont(a).square()), fe(a).square().value);
        }

        #[test]
        fn mont_neg_matches_default(a in 0u64..17) {
            prop_assert_eq!(canonical(-fe_mont(a)), (-fe(a)).value);
        }

        #[test]
        fn mont_inverse_matches_default(a in 1u64..17) {
            prop_assert_eq!(canonical(fe_mont(a).inverse()), fe(a).inverse().value);
        }

        #[test]
        fn mont_div_matches_default(a in 0u64..17, b in 1u64..17) {
            prop_assert_eq!(canonical(fe_mont(a) / fe_mont(b)), (fe(a) / fe(b)).value);
        }
    }

    #[test]
    fn pow_zero_exponent_is_one() {
        assert_eq!(fe(5).pow(U256::ZERO).value, U256::from(1));
    }

    #[test]
    fn pow_one_exponent_is_self() {
        assert_eq!(fe(5).pow(U256::ONE).value, fe(5).value);
    }

    proptest! {
        #[test]
        fn pow_matches_repeated_mul(a in 0u64..17, exp in 0u32..8) {
            let mut expected = fe(1);
            for _ in 0..exp {
                expected = expected * fe(a);
            }

            prop_assert_eq!(fe(a).pow(U256::from(exp as u64)).value, expected.value);
        }

        #[test]
        fn pow_result_is_canonical(a in 0u64..17, exp in 0u32..8) {
            prop_assert!(fe(a).pow(U256::from(exp as u64)).value < Mod17::MODULUS);
        }

        #[test]
        fn mont_pow_matches_default(a in 0u64..17, exp in 0u32..8) {
            prop_assert_eq!(
                canonical(fe_mont(a).pow(U256::from(exp as u64))),
                fe(a).pow(U256::from(exp as u64)).value
            );
        }

        #[test]
        fn fermats_little_theorem(a in 1u64..17) {
            // a^(p-1) == 1 (mod p) for prime p and a not divisible by p.
            prop_assert_eq!(fe(a).pow(U256::from(16)).value, U256::from(1));
        }
    }

}
