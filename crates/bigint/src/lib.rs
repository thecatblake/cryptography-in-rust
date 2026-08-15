use std::ops::{Index, Add, Sub, Mul, Div, Rem, AddAssign, SubAssign, Shl, ShlAssign, Shr, ShrAssign};

use std::cmp::Ordering;

use field_core::{EuclideanRepr, FpRepr, WideInt};

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

impl<const N: usize> FpRepr for Uint<N> {
    const ZERO: Self = Self::ZERO;
    const BITS: usize = N * 64;

    fn bit(&self, i: usize) -> bool {
        Uint::bit(self, i)
    }
}

// Sub/wrapping_mul both wrap mod 2^(64*N) (two's-complement semantics),
// so is_negative's top-bit check is the same trick used to recover a
// canonical positive residue from a "negative" intermediate Bezout
// coefficient in field_core::gcd_inverse.
impl<const N: usize> EuclideanRepr for Uint<N> {
    const ONE: Self = Self::ONE;

    fn wrapping_sub(self, rhs: Self) -> Self {
        self - rhs
    }

    fn wrapping_mul(self, rhs: Self) -> Self {
        Uint::wrapping_mul(self, rhs)
    }

    fn div_rem(self, rhs: Self) -> (Self, Self) {
        Uint::div_rem(self, rhs)
    }

    fn is_negative(self) -> bool {
        self.bit(N * 64 - 1)
    }
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
    // a*a is symmetric across the limb grid's diagonal (a[i]*a[j] ==
    // a[j]*a[i]), so each off-diagonal product only needs to be computed
    // once and doubled, instead of twice as full Mul would. That's 10
    // single-limb multiplications (4 diagonal + 6 cross) instead of 16.
    pub fn square(self) -> U512 {
        let a = self;

        // Cross sum: sum of a[i]*a[j] for i < j, positioned at limb i+j.
        let mut cross = U512::ZERO;
        for i in 0..4 {
            let mut row = [0u64; 8];
            let mut carry = 0u64;

            for j in (i + 1)..4 {
                let sum = a[i] as u128 * a[j] as u128 + carry as u128;

                row[i + j] = sum as u64;
                carry = (sum >> 64) as u64;
            }
            row[i + 4] = carry;

            cross = cross + U512 { limbs: row };
        }

        // a^2 = diag + 2*cross, and a^2 < 2^512, so 2*cross < 2^512 and
        // this left shift can't lose bits off the top.
        let cross_doubled = cross << 1;

        // Diagonal terms a[i]^2 land at limb pair [2i, 2i+1]; consecutive
        // pairs never overlap, so no carry propagation is needed between them.
        let mut diag = [0u64; 8];
        for i in 0..4 {
            let sq = a[i] as u128 * a[i] as u128;
            diag[2 * i] = sq as u64;
            diag[2 * i + 1] = (sq >> 64) as u64;
        }

        cross_doubled + U512 { limbs: diag }
    }

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
            base = base.square().resize();
            e >>= 1;
        }

        result
    }
}

// Only U256 gets this (not a blanket impl over Uint<N>): stable Rust has no
// const-generic arithmetic to spell "the double-width Uint<2*N>" generically,
// so each doubling pair is written out by hand -- same reasoning as U256's
// own Mul/square impls above, and the same pattern field-core's
// impl_wide_int! macro follows for primitive pairs (u64 => u128, etc). This
// is what lets a WideFieldConfig (field-core) be implemented directly on
// U256, e.g. for a 256-bit prime with a fast reduction trick: add/sub/neg
// come out of WideEuclideanBackend for free via the widen()/narrow() round
// trip through U512, and the implementor only has to write `mul`.
impl WideInt for U256 {
    type Wide = U512;

    fn widen(self) -> U512 {
        self.resize()
    }

    fn narrow(wide: U512) -> U256 {
        wide.resize()
    }

    fn wide_mul(self, other: U256) -> U512 {
        self * other
    }

    fn from_u8(v: u8) -> U256 {
        U256::from_u64(v as u64)
    }
}

impl<const N: usize> Shl<usize> for Uint<N> {
    type Output = Self;

    // Shifts spanning more than one limb (rhs >= 64) need to move whole
    // limbs first, then apply the remaining sub-64 bit shift across each
    // adjacent limb pair. div_rem relies on this: it aligns divisor with
    // dividend by shifting by (dividend.bits() - divisor.bits()), which
    // routinely exceeds 63 once dividend and divisor differ by more than
    // one limb in magnitude (e.g. a U512 product reduced by a U256 modulus).
    fn shl(self, rhs: usize) -> Self::Output {
        assert!(rhs < N * 64);

        if rhs == 0 {
            return self;
        }

        let limb_shift = rhs / 64;
        let bit_shift = rhs % 64;

        let mut result = [0u64; N];

        for i in (limb_shift..N).rev() {
            let mut v = self[i - limb_shift] << bit_shift;

            if bit_shift != 0 && i > limb_shift {
                v |= self[i - limb_shift - 1] >> (64 - bit_shift);
            }

            result[i] = v;
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

    // See Shl's comment: rhs can span multiple limbs, so shift whole limbs
    // first and then the remaining sub-64 bit shift per adjacent limb pair.
    fn shr(self, rhs: usize) -> Self::Output {
        assert!(rhs < N * 64);

        if rhs == 0 {
            return self;
        }

        let limb_shift = rhs / 64;
        let bit_shift = rhs % 64;

        let mut result = [0u64; N];

        for i in 0..(N - limb_shift) {
            let src = i + limb_shift;
            let mut v = self[src] >> bit_shift;

            if bit_shift != 0 && src + 1 < N {
                v |= self[src + 1] << (64 - bit_shift);
            }

            result[i] = v;
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

    // Schoolbook multiply, truncated to N limbs (i.e. mod 2^(64*N)) instead
    // of widening -- used by EuclideanRepr's extended-GCD Bezout
    // coefficients, which are only ever consumed mod 2^(64*N) and are
    // allowed to wrap the same way Add/Sub above do.
    pub fn wrapping_mul(self, rhs: Self) -> Self {
        let mut result = [0u64; N];

        for i in 0..N {
            if self[i] == 0 {
                continue;
            }

            let mut carry = 0u128;
            for j in 0..(N - i) {
                let sum = self[i] as u128 * rhs[j] as u128 + result[i + j] as u128 + carry;
                result[i + j] = sum as u64;
                carry = sum >> 64;
            }
        }

        Uint { limbs: result }
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

    #[test]
    fn square_zero() {
        assert_eq!(U256::ZERO.square(), U512::ZERO);
    }

    #[test]
    fn square_one() {
        assert_eq!(U256::ONE.square(), U512::ONE);
    }

    #[test]
    fn square_propagates_carries_across_limbs() {
        // (2^256 - 1)^2 exercises carries through both the cross-sum
        // accumulation and the final cross+diag addition.
        let all_ones = U256::from_limbs([u64::MAX; 4]);

        assert_eq!(all_ones.square(), all_ones * all_ones);
    }

    proptest! {
        #[test]
        fn square_matches_mul(a in any::<u64>()) {
            let ua = U256::from(a);

            prop_assert_eq!(ua.square(), ua * ua);
        }

        #[test]
        fn square_matches_u64(a in any::<u64>()) {
            let ua = U256::from(a);
            let expected = (a as u128) * (a as u128);

            prop_assert_eq!(ua.square().low_u128(), expected);
        }
    }

    #[test]
    fn shl_within_single_limb_matches_old_behavior() {
        assert_eq!(U256::from(1u64) << 3, U256::from(8u64));
    }

    #[test]
    fn shl_by_zero_is_identity() {
        let a = U256::from_limbs([1, 2, 3, 4]);
        assert_eq!(a << 0, a);
    }

    #[test]
    fn shl_across_one_limb_boundary() {
        // 1 << 64 must land exactly on limb 1, not panic like the old
        // single-limb-only Shl (assert!(rhs < 64)) would have.
        assert_eq!(U256::from(1u64) << 64, U256::from_limbs([0, 1, 0, 0]));
    }

    #[test]
    fn shl_across_multiple_limbs_with_sub_limb_remainder() {
        // 1 << 130 = limb_shift 2, bit_shift 2 -> bit 130 set, i.e. limb 2 = 0b100.
        assert_eq!(U256::from(1u64) << 130, U256::from_limbs([0, 0, 0b100, 0]));
    }

    #[test]
    fn shl_carries_bits_across_the_shifted_boundary() {
        // Bit 63 shifted by 65 lands at bit 128 (limb 2, offset 0): the
        // in-limb shift wraps a bit up past the immediately-next limb.
        let a = U256::from_limbs([1u64 << 63, 0, 0, 0]);
        assert_eq!(a << 65, U256::from_limbs([0, 0, 1, 0]));
    }

    #[test]
    fn shr_across_multiple_limbs_with_sub_limb_remainder() {
        let a = U256::from_limbs([0, 0, 0b100, 0]);
        assert_eq!(a >> 130, U256::from(1u64));
    }

    #[test]
    fn shr_carries_bits_across_the_shifted_boundary() {
        // Bit 65 shifted right by 65 lands at bit 0.
        let a = U256::from_limbs([0, 2, 0, 0]);
        assert_eq!(a >> 65, U256::from(1u64));
    }

    #[test]
    fn shl_shr_roundtrip_across_limb_boundary() {
        // Kept under 64 bits so shifting left by 130 (194 bits total) can't
        // lose bits off the top of a U256, letting the round trip hold.
        let a = U256::from_limbs([0x1234_5678_9abc_def0, 0, 0, 0]);
        assert_eq!((a << 130) >> 130, a);
    }

    proptest! {
        #[test]
        fn shl_matches_u128_for_small_shifts(a in any::<u64>(), shift in 0usize..128) {
            let ua = U256::from(a);
            let expected = ((a as u128) << shift) as u64;

            prop_assert_eq!((ua << shift).low_u64(), expected);
        }

        #[test]
        fn shr_matches_u128_for_small_shifts(a in any::<u64>(), shift in 0usize..128) {
            let ua = U256::from(a);
            // a fits in the low limb, so shifting right by >= 64 always yields 0.
            let expected = if shift >= 64 { 0 } else { a >> shift };

            prop_assert_eq!((ua >> shift).low_u64(), expected);
        }
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
    fn div_rem_wide_bit_length_gap_regression() {
        // dividend.bits() - divisor.bits() here is ~193, forcing div_rem's
        // `divisor <<= shift` past a single limb. The old Shl (assert!(rhs
        // < 64)) panicked on this; this is exactly the shape DefaultBackend's
        // `product % modulus` produces for a real-sized (non-tiny) prime.
        let dividend = U512::from_limbs([0, 0, 0, 1, 0, 0, 0, 0]); // 2^192
        let divisor: U512 = U256::from(3u64).resize();

        let (q, r) = dividend.div_rem(divisor);

        // 2 == -1 (mod 3), so 2^192 == (-1)^192 == 1 (mod 3).
        assert_eq!(r, U512::from(1u64));
        assert!(r < divisor);
        // U512 has no Mul, so verify q*3 + r == dividend via repeated addition.
        assert_eq!(q + q + q + r, dividend);
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
