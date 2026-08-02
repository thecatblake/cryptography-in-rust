use bigint::{U256, U512};
use bigint::math::extended_gcd;
use std::ops::{Add, Sub, Mul, Div, Neg};
use std::marker::PhantomData;

// FpConfig associates a modulus with the type rather than each value.
// This makes field elements with different moduli different Rust types,
// and we can impose the operation under the same modulus at the compile time

pub trait FpConfig {
    const MODULUS: U256;
}

pub struct Fp<C: FpConfig> {
    value: bigint::U256,
    _marker: PhantomData<C>
}

impl<C: FpConfig> Fp<C> {
    pub fn new(value: bigint::U256) -> Self {
        Fp { value, _marker: PhantomData }
    }
}

impl<C: FpConfig> Add for Fp<C> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut sum = self.value + rhs.value;

        if sum >= C::MODULUS {
            sum -= C::MODULUS;
        }

        Fp::new(sum)
    }
}

impl<C: FpConfig> Sub for Fp<C> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let diff = if self.value >= rhs.value {
            self.value - rhs.value
        } else {
            (self.value + C::MODULUS) - rhs.value
        };

        Fp::new(diff)
    }
}

impl<C: FpConfig> Mul for Fp<C> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let product: U512 = self.value * rhs.value;
        let modulus: U512 = C::MODULUS.resize();

        Fp::new((product % modulus).resize())
    }
}

impl<C: FpConfig> Fp<C> {
    // Multiplicative inverse mod C::MODULUS, via the extended Euclidean algorithm.
    // ax + bp == ax == 1 mod p where p is prime
    pub fn inverse(self) -> Self {
        assert!(self.value != U256::ZERO, "cannot invert zero in a field");

        let egcd = extended_gcd(self.value, C::MODULUS);

        let inv = if egcd.x_neg {
            C::MODULUS - egcd.x
        } else {
            egcd.x
        };

        Fp::new(inv)
    }
}

impl<C: FpConfig> Div for Fp<C> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.inverse()
    }
}

impl<C: FpConfig> Neg for Fp<C> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let value = if self.value == U256::ZERO {
            U256::ZERO
        } else {
            C::MODULUS - self.value
        };

        Fp::new(value)
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

    type F17 = Fp<Mod17>;

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
}
