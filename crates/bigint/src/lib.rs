use std::ops::{Index, Add, Sub, Mul};

#[derive(Clone, Debug, Eq, PartialEq, Copy)]
pub struct Uint<const N: usize>{
    limbs: [u64; N],
}

impl<const N: usize> Uint<N> {
    pub const fn low_u64(&self) -> u64 {
        self.limbs[0]
    }

    pub fn low_u128(&self) -> u128 {
        ((self[1] as u128) << 64) | (self[0] as u128)
    }
}

pub type U256 = Uint<4>;
pub type U512 = Uint<8>;

impl<const N: usize> Uint<N> {
    pub const ZERO: Self = Self { limbs: [0; N] };
}

impl<const N: usize> From<u64> for Uint<N> {
    fn from(value: u64) -> Self {
        let mut limbs = [0u64; N];

        if N > 0 {
            limbs[0] = value;
        }

        Self { limbs }
    }
}

impl<const N: usize> Index<usize> for Uint<N> {
    type Output = u64;

    fn index(&self, index: usize) -> &Self::Output {
        &self.limbs[index]
    }
}

impl<const N:usize> Add for Uint<N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut result = [0u64; N];
        let mut carry = 0u128;

        for i in 0..N {
            let sum = self[i] as u128 
                + rhs[i] as u128
                + carry;
            result[i] = sum as u64;
            carry = sum >> 64;
        }

        Uint { limbs: result }
    }
}

impl<const N:usize> Sub for Uint<N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut result = [0u64; N];
        let mut borrow = 0u64;

        for i in 0..N {
            let (tmp, o1) = self[i].overflowing_sub(rhs[i]);
            let (res, o2) = tmp.overflowing_sub(borrow);
            borrow = (o1 || o2) as u64;
            result[i] = res;
        }

        Uint { limbs: result }
    }
}

impl Mul for U256 {
    type Output = U512;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = [0u64; 8];

        for i in 0..4 {
            for j in 0..4 {
                let mul = self[i] as u128
                    * rhs[j] as u128;

                result[i + j] += mul as u64;
                result[i + j + 1] += (mul >> 64) as u64;
            }
        }

        U512 { limbs: result }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn add_zero(a in any::<u64>()) {
            let a = U256::from(a);
            prop_assert_eq!(a + U256::ZERO, a);
        }

        #[test]
        fn add_commutative(a in any::<u64>(), b in any::<u64>()) {
            let a = U256::from(a);
            let b = U256::from(b);

            prop_assert_eq!(a + b, b + a);
        }

        #[test]
        fn add_associative(
            a in any::<u64>(),
            b in any::<u64>(),
            c in any::<u64>(),
        ) {
            let a = U256::from(a);
            let b = U256::from(b);
            let c = U256::from(c);

            prop_assert_eq!((a + b) + c, a + (b + c));
        }

        #[test]
        fn sub_self(a in any::<u64>()) {
            let a = U256::from(a);

            prop_assert_eq!(a - a, U256::ZERO);
        }

        #[test]
        fn subtraction_roundtrip(
            a in any::<u64>(),
            b in any::<u64>(),
        ) {
            prop_assume!(a >= b);

            let a = U256::from(a);
            let b = U256::from(b);

            prop_assert_eq!((a - b) + b, a);
        }
    
        #[test]
        fn sub_zero(a in any::<u64>()) {
            let a = U256::from(a);

            prop_assert_eq!(a - U256::ZERO, a);
        }
        #[test]
        fn add_matches_u64(a in any::<u64>(), b in any::<u64>()) {
            let ua = U256::from(a);
            let ub = U256::from(b);

            let expected = a.wrapping_add(b);

            prop_assert_eq!(
                (ua + ub).low_u64(),
                expected
            );
        }

        #[test]
        fn sub_matches_u64(a in any::<u64>(), b in any::<u64>()) {
            let ua = U256::from(a);
            let ub = U256::from(b);

            let expected = a.wrapping_sub(b);

            prop_assert_eq!(
                (ua - ub).low_u64(),
                expected
            );
        }

        #[test]
        fn mul_unit(a in any::<u64>()) {
            let ua = U256::from(a);
            let unit = U256::from(1);

            prop_assert_eq!(
                ua * unit,
                U512::from(a)
            );
        }

        #[test]
        fn mul_zero(a in any::<u64>()) {
            let ua = U256::from(a);
            let zero = U256::from(0);

            prop_assert_eq!(
                ua * zero,
                U512::from(0)
            );
        }

        #[test]
        fn mul_commutative(a in any::<u64>(), b in any::<u64>()) {
            let a = U256::from(a);
            let b = U256::from(b);

            prop_assert_eq!(a * b, b * a);
        }

        #[test]
        fn mul_matches_u64(a in any::<u64>(), b in any::<u64>()) {
            let ua = U256::from(a);
            let ub = U256::from(b);

            let expected = (a as u128) * (b as u128);

            let result = ua * ub;

            assert_eq!(result.low_u128(), expected);
        }
    }
}
