use crate::Uint;

pub fn gcd<const N: usize>(a: Uint<N>, b: Uint<N>) -> Uint<N> {
    if a < b  {
        gcd(b, a)
    }
    else if b == Uint::ZERO {
        a
    }
    else if a == Uint::ZERO {
        b
    }
    else {
        let r = a % b;
        gcd(b, r)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ExtendedGcd<const N: usize> {
    pub gcd: Uint<N>,
    pub x: Uint<N>,
    pub x_neg: bool,
    pub y: Uint<N>,
    pub y_neg: bool,
}

fn mul_wrapping<const N: usize>(a: Uint<N>, b: Uint<N>) -> Uint<N> {
    let mut result = [0u64; N];

    for i in 0..N {
        if a[i] == 0 {
            continue;
        }

        let mut carry = 0u128;
        for j in 0..(N - i) {
            let sum = a[i] as u128 * b[j] as u128
                + result[i + j] as u128
                + carry;
            result[i + j] = sum as u64;
            carry = sum >> 64;
        }
    }

    Uint { limbs: result }
}

fn is_negative<const N: usize>(v: Uint<N>) -> bool {
    v.bit(N * 64 - 1)
}

fn negate<const N: usize>(v: Uint<N>) -> Uint<N> {
    Uint::ZERO - v
}

pub fn extended_gcd<const N: usize>(a: Uint<N>, b: Uint<N>) -> ExtendedGcd<N> {
    let mut old_r = a;
    let mut r = b;
    let mut old_s = Uint::ONE;
    let mut s = Uint::ZERO;
    let mut old_t = Uint::ZERO;
    let mut t = Uint::ONE;

    while r != Uint::ZERO {
        let (q, rem) = old_r.div_rem(r);

        old_r = r;
        r = rem;

        let new_s = old_s - mul_wrapping(q, s);
        old_s = s;
        s = new_s;

        let new_t = old_t - mul_wrapping(q, t);
        old_t = t;
        t = new_t;
    }

    let (x, x_neg) = if is_negative(old_s) {
        (negate(old_s), true)
    } else {
        (old_s, false)
    };

    let (y, y_neg) = if is_negative(old_t) {
        (negate(old_t), true)
    } else {
        (old_t, false)
    };

    ExtendedGcd { gcd: old_r, x, x_neg, y, y_neg }
}

// Multiplicative inverse mod C::MODULUS, via the extended Euclidean algorithm.
// ax + bp == ax == 1 mod p where p is prime
pub fn mod_inv<const N: usize>(a: Uint<N>, n: Uint<N>) -> Uint<N> {
    let ext = extended_gcd(a, n);
    if ext.x_neg {
        n - ext.x
    } else {
        ext.x
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::U256;
    use proptest::prelude::*;

    fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
        while b != 0 {
            let t = b;
            b = a % b;
            a = t;
        }
        a
    }

    #[test]
    fn gcd_zero_zero() {
        assert_eq!(gcd(U256::ZERO, U256::ZERO), U256::ZERO);
    }

    #[test]
    fn gcd_a_zero() {
        let a = U256::from(42);
        assert_eq!(gcd(a, U256::ZERO), a);
    }

    #[test]
    fn gcd_zero_b() {
        let b = U256::from(42);
        assert_eq!(gcd(U256::ZERO, b), b);
    }

    fn extended_gcd_i128(a: i128, b: i128) -> (i128, i128, i128) {
        let (mut old_r, mut r) = (a, b);
        let (mut old_s, mut s) = (1i128, 0i128);
        let (mut old_t, mut t) = (0i128, 1i128);

        while r != 0 {
            let q = old_r / r;

            let new_r = old_r - q * r;
            old_r = r;
            r = new_r;

            let new_s = old_s - q * s;
            old_s = s;
            s = new_s;

            let new_t = old_t - q * t;
            old_t = t;
            t = new_t;
        }

        (old_r, old_s, old_t)
    }

    fn signed(mag: U256, neg: bool) -> i128 {
        let mag = mag.low_u64() as i128;
        if neg { -mag } else { mag }
    }

    #[test]
    fn extended_gcd_zero_zero() {
        let r = extended_gcd(U256::ZERO, U256::ZERO);
        assert_eq!(r.gcd, U256::ZERO);
        assert_eq!(r.x, U256::ONE);
        assert!(!r.x_neg);
        assert_eq!(r.y, U256::ZERO);
    }

    #[test]
    fn extended_gcd_a_zero() {
        let a = U256::from(42);
        let r = extended_gcd(a, U256::ZERO);
        assert_eq!(r.gcd, a);
        assert_eq!(r.x, U256::ONE);
        assert!(!r.x_neg);
        assert_eq!(r.y, U256::ZERO);
    }

    #[test]
    fn extended_gcd_zero_b() {
        let b = U256::from(42);
        let r = extended_gcd(U256::ZERO, b);
        assert_eq!(r.gcd, b);
        assert_eq!(r.x, U256::ZERO);
        assert_eq!(r.y, U256::ONE);
        assert!(!r.y_neg);
    }

    #[test]
    fn extended_gcd_12_8() {
        let r = extended_gcd(U256::from(12), U256::from(8));
        assert_eq!(r.gcd, U256::from(4));
        assert_eq!(r.x, U256::from(1));
        assert!(!r.x_neg);
        assert_eq!(r.y, U256::from(1));
        assert!(r.y_neg);
    }

    #[test]
    fn mod_inv_3_11() {
        assert_eq!(mod_inv(U256::from(3), U256::from(11)), U256::from(4));
    }

    #[test]
    fn mod_inv_one() {
        assert_eq!(mod_inv(U256::from(1), U256::from(7)), U256::from(1));
    }

    #[test]
    fn mod_inv_rsa_example() {
        // e = 17, phi = 3120 -> d = 2753 (17 * 2753 = 46801 = 15*3120 + 1)
        assert_eq!(mod_inv(U256::from(17), U256::from(3120)), U256::from(2753));
    }

    #[test]
    fn mod_inv_reduces_a_greater_than_n() {
        // 40 mod 7 = 5, and 5 * 3 = 15 = 1 mod 7
        assert_eq!(mod_inv(U256::from(40), U256::from(7)), U256::from(3));
    }

    proptest! {
        #[test]
        fn gcd_commutative(a in any::<u64>(), b in any::<u64>()) {
            let ua = U256::from(a);
            let ub = U256::from(b);

            prop_assert_eq!(gcd(ua, ub), gcd(ub, ua));
        }

        #[test]
        fn gcd_matches_u64(a in any::<u64>(), b in any::<u64>()) {
            let ua = U256::from(a);
            let ub = U256::from(b);

            prop_assert_eq!(gcd(ua, ub).low_u64(), gcd_u64(a, b));
        }

        #[test]
        fn gcd_divides_both(a in 1u64.., b in 1u64..) {
            let ua = U256::from(a);
            let ub = U256::from(b);

            let g = gcd(ua, ub);

            prop_assert_eq!(ua % g, U256::ZERO);
            prop_assert_eq!(ub % g, U256::ZERO);
        }

        #[test]
        fn extended_gcd_matches_i128(a in any::<u64>(), b in any::<u64>()) {
            let r = extended_gcd(U256::from(a), U256::from(b));

            let (expected_gcd, expected_x, expected_y) =
                extended_gcd_i128(a as i128, b as i128);

            prop_assert_eq!(r.gcd.low_u64() as i128, expected_gcd);
            prop_assert_eq!(signed(r.x, r.x_neg), expected_x);
            prop_assert_eq!(signed(r.y, r.y_neg), expected_y);
        }

        #[test]
        fn extended_gcd_bezout_identity(a in 1u64.., b in 1u64..) {
            let ua = U256::from(a);
            let ub = U256::from(b);

            let r = extended_gcd(ua, ub);

            let ax = signed(ua, false) as i128 * signed(r.x, r.x_neg);
            let by = signed(ub, false) as i128 * signed(r.y, r.y_neg);

            prop_assert_eq!(ax + by, r.gcd.low_u64() as i128);
        }

        #[test]
        fn mod_inv_is_multiplicative_inverse(a in 1u64.., n in 2u64..) {
            prop_assume!(gcd_u64(a, n) == 1);

            let inv = mod_inv(U256::from(a), U256::from(n));

            prop_assert!(inv < U256::from(n));
            prop_assert_eq!((a as u128 * inv.low_u64() as u128) % (n as u128), 1);
        }
    }
}