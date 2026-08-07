use std::ops::{Index, Add, Sub, Mul, Div, Rem, AddAssign, SubAssign, Shl, ShlAssign, Shr, ShrAssign};

use std::cmp::Ordering;

pub mod math;

#[derive(Clone, Debug, Eq, PartialEq, Copy)]
pub struct Uint<const N: usize>{
    limbs: [u64; N],
}

pub type U256 = Uint<4>;
pub type U512 = Uint<8>;

impl<const N: usize> Uint<N> {
    pub const ONE: Self = {
        let mut limbs = [0; N];
        limbs[0] = 1;
        Self { limbs }
    };

    pub const fn low_u64(&self) -> u64 {
        self.limbs[0]
    }

    pub const fn from_u64(value: u64) -> Self {
        let mut limbs = [0u64; N];

        if N > 0 {
            limbs[0] = value;
        }

        Self { limbs }
    }

    pub const fn from_limbs(limbs: [u64; N]) -> Self {
        Self { limbs }
    }

    pub fn low_u128(&self) -> u128 {
        ((self[1] as u128) << 64) | (self[0] as u128)
    }

    // Zero-extends into a wider width, or truncates the high limbs when narrowing.
    pub fn resize<const M: usize>(&self) -> Uint<M> {
        let mut limbs = [0u64; M];
        let len = N.min(M);

        limbs[..len].copy_from_slice(&self.limbs[..len]);

        Uint { limbs }
    }

    pub fn bit(&self, i: usize) -> bool {
        let limb = i / 64;
        let offset = i % 64;

        ((self[limb] >> offset) & 1) == 1
    }

    pub fn set_bit(&mut self, i: usize) {
        let limb = i / 64;
        let offset = i % 64;

        self.limbs[limb] |= 1 << offset;
    }

    pub fn bits(&self) -> usize {
        for limb in (0..N).rev() {
            let x = self.limbs[limb];
            if x != 0 {
                return limb * 64 + (63 - x.leading_zeros() as usize);
            }
        }
        0
    }
}


impl<const N: usize> Ord for Uint<N> {
    fn cmp(&self, other: &Self) -> Ordering {
        for i in (0..N).rev() {
            match self.limbs[i].cmp(&other.limbs[i]) {
                Ordering::Equal => continue,
                ord => return ord,
            }
        }

        Ordering::Equal
    }
}

impl<const N: usize> PartialOrd for Uint<N> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<const N: usize> Uint<N> {
    pub const ZERO: Self = Self { limbs: [0; N] };
}

impl<const N: usize> From<u64> for Uint<N> {
    fn from(value: u64) -> Self {
        Self::from_u64(value)
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

impl<const N: usize> AddAssign for Uint<N> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
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

impl<const N: usize> SubAssign for Uint<N> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Mul for U256 {
    type Output = U512;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = U512::ZERO;

        for i in 0..4 {
            let mut row = [0u64; 8];
            let mut carry = 0u64;

            for j in 0..4 {
                let sum = self[i] as u128 * rhs[j] as u128 + carry as u128;

                row[i + j] = sum as u64;
                carry = (sum >> 64) as u64;
            }
            row[i + 4] = carry;

            result = result + U512 { limbs: row };
        }

        result
    }
}

impl U256 {
    // Square-and-multiply, wrapping mod 2^256 (matches the widening Mul
    // above, truncated back down via resize on each step).
    pub fn pow(self, exp: Self) -> Self {
        let mut result = Self::ONE;
        let mut base = self;
        let mut e = exp;

        while e != Self::ZERO {
            if e.bit(0) {
                result = (result * base).resize();
            }
            base = (base * base).resize();
            e >>= 1;
        }

        result
    }
}

impl<const N: usize> Shl<usize> for Uint<N> {
    type Output = Self;

    fn shl(self, rhs: usize) -> Self::Output {
        assert!(rhs < 64);

        let mut result = [0u64; N];
        let mut carry = 0u64;

        for i in 0..N {
            let new_carry = if rhs == 0 {
                0
            } else {
                self[i] >> (64 - rhs)
            };

            result[i] = (self[i] << rhs) | carry;
            carry = new_carry;
        }

        Uint { limbs: result }
    }
}

impl<const N: usize> ShlAssign<usize> for Uint<N> {
    fn shl_assign(&mut self, rhs: usize) {
        *self = *self << rhs;
    }
}

impl<const N: usize> Shr<usize> for Uint<N> {
    type Output = Self;

    fn shr(self, rhs: usize) -> Self::Output {
        assert!(rhs < 64);

        let mut result = [0u64; N];
        let mut carry = 0u64;

        for i in (0..N).rev() {
            let new_carry = if rhs == 0 {
                0
            } else {
                self[i] << (64 - rhs)
            };

            result[i] = (self[i] >> rhs) | carry;
            carry = new_carry;
        }

        Uint { limbs: result }
    }
}

impl<const N: usize> ShrAssign<usize> for Uint<N> {
    fn shr_assign(&mut self, rhs: usize) {
        *self = *self >> rhs;
    }
}

impl<const N: usize> Uint<N> {
    pub fn div_rem(self, rhs: Self) -> (Self, Self) {
        assert!(rhs != Self::ZERO);

        if self < rhs {
            return (Self::ZERO, self);
        }

        if self == rhs {
            return (Self::ONE, Self::ZERO);
        }

        let mut dividend = self;
        let mut divisor = rhs;

        let mut quotient = Self::ZERO;

        let shift = dividend.bits() - divisor.bits();

        divisor <<= shift;

        for i in (0..=shift).rev() {
            if dividend >= divisor {
                dividend -= divisor;
                quotient.set_bit(i);
            }

            divisor >>= 1;
        }

        (quotient, dividend)
    }
}

impl<const N: usize> Div for Uint<N> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self.div_rem(rhs).0
    }
}

impl<const N: usize> Rem for Uint<N> {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        self.div_rem(rhs).1
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

    #[test]
    fn mul_propagates_carries_across_limbs() {
        // (2^256 - 1) * 2 = 2^257 - 2, which requires a carry to ripple
        // through every limb of the low half and into the high half.
        // The old accumulate-without-carry Mul impl overflowed a limb here.
        let all_ones = U256::from_limbs([u64::MAX; 4]);
        let two = U256::from_u64(2);

        let expected = U512::from_limbs([
            u64::MAX - 1,
            u64::MAX,
            u64::MAX,
            u64::MAX,
            1,
            0,
            0,
            0,
        ]);

        assert_eq!(all_ones * two, expected);
    }

    proptest! {
        #[test]
        fn div_by_one(a in any::<u64>()) {
            let ua = U256::from(a);
            let one = U256::from(1);

            prop_assert_eq!(ua / one, ua);
        }

        #[test]
        fn div_self(a in any::<u64>()) {
            prop_assume!(a != 0);

            let ua = U256::from(a);

            prop_assert_eq!(ua / ua, U256::ONE);
        }

        #[test]
        fn zero_div(a in any::<u64>()) {
            prop_assume!(a != 0);

            let ua = U256::from(a);

            prop_assert_eq!(U256::ZERO / ua, U256::ZERO);
        }

        #[test]
        fn div_matches_u64(a in any::<u64>(), b in any::<u64>()) {
            prop_assume!(b != 0);

            let ua = U256::from(a);
            let ub = U256::from(b);

            let expected = a / b;

            prop_assert_eq!(
                (ua / ub).low_u64(),
                expected
            );
        }

        #[test]
        fn rem_matches_u64(a in any::<u64>(), b in any::<u64>()) {
            prop_assume!(b != 0);

            let ua = U256::from(a);
            let ub = U256::from(b);

            let expected = a % b;

            prop_assert_eq!(
                (ua % ub).low_u64(),
                expected
            );
        }

        #[test]
        fn rem_zero(a in any::<u64>()) {
            prop_assume!(a != 0);

            let ua = U256::from(a);

            prop_assert_eq!(ua % ua, U256::ZERO);
        }

        #[test]
        fn div_rem_roundtrip(a in any::<u64>(), b in any::<u64>()) {
            prop_assume!(b != 0);

            let ua = U256::from(a);
            let ub = U256::from(b);

            let (q, r) = ua.div_rem(ub);

            prop_assert_eq!(q * ub + U512::from(r.low_u64()), U512::from(a));
            prop_assert!(r < ub);
        }
    }

    #[test]
    fn pow_zero_exponent_is_one() {
        assert_eq!(U256::from(5).pow(U256::ZERO), U256::ONE);
    }

    #[test]
    fn pow_zero_base_is_zero() {
        assert_eq!(U256::ZERO.pow(U256::from(5)), U256::ZERO);
    }

    #[test]
    fn pow_one_exponent_is_self() {
        let a = U256::from(12345);
        assert_eq!(a.pow(U256::ONE), a);
    }

    #[test]
    fn pow_wraps_mod_2_256() {
        // 2^256 == 0 (mod 2^256), matching wrapping semantics.
        assert_eq!(U256::from(2).pow(U256::from(256)), U256::ZERO);
    }

    proptest! {
        #[test]
        fn pow_matches_u64_when_no_overflow(base in 0u64..1000, exp in 0u32..5) {
            let ua = U256::from(base);
            let expected = base.pow(exp);

            prop_assert_eq!(ua.pow(U256::from(exp as u64)).low_u64(), expected);
        }

        #[test]
        fn pow_matches_repeated_mul(base in any::<u64>(), exp in 0u32..12) {
            let ua = U256::from(base);

            let mut expected = U256::ONE;
            for _ in 0..exp {
                expected = (expected * ua).resize();
            }

            prop_assert_eq!(ua.pow(U256::from(exp as u64)), expected);
        }
    }
}
