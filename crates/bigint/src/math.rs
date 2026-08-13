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
    }
}
