use field_core::{Fp2, Fp2Config, Fp6, Fp6Config, Fp12, Fp12Config};

use crate::{Mersenne31, Mersenne31Config};

// Fp2/Fp6/Fp12 tower over Mersenne31. Unlike babybear's Montgomery backend,
// Mersenne31Backend represents values in plain canonical form, so there's no
// to_mont/pow_mont_* compile-time machinery to lean on here (that machinery
// is specific to MontFieldConfig) -- every constant below is just the
// canonical residue found by the accompanying offline derivation script
// (small-non-residue search + modular exponentiation over the exact same
// QuadExt/CubicExt tower arithmetic this crate builds on), the same
// "verified offline" convention field-core's own Fp12 fixtures use for
// GAMMA/FROBENIUS_COEFF (no pow_mont_fp6 helper exists at any level here).
impl Fp2Config for Mersenne31Config {
    type Base = Mersenne31;

    // 3 is the smallest quadratic non-residue mod p = 2^31-1.
    const BETA: Self::Base = Mersenne31::new(3);

    // BETA^((p-1)/2) mod p.
    const FROBENIUS_COEFF: Self::Base = Mersenne31::new(2_147_483_646);
}

impl Fp6Config for Mersenne31Config {
    // 1 + 2u is the smallest cubic non-residue in Fp2 found by the search.
    const XI: Fp2<Self> = Fp2 { c0: Mersenne31::new(1), c1: Mersenne31::new(2) };

    // XI^((p-1)/3) and XI^(2(p-1)/3) in Fp2.
    const FROBENIUS_COEFF_C1: Fp2<Self> =
        Fp2 { c0: Mersenne31::new(2_044_269_142), c1: Mersenne31::new(697_200_390) };
    const FROBENIUS_COEFF_C2: Fp2<Self> =
        Fp2 { c0: Mersenne31::new(1_409_171_703), c1: Mersenne31::new(647_605_448) };
}

impl Fp12Config for Mersenne31Config {
    // v (Fp6's own generator, embedded via c1 with c0 = c2 = 0) is the
    // smallest quadratic non-residue in Fp6 found by the search.
    const GAMMA: Fp6<Self> = Fp6 {
        c0: Fp2 { c0: Mersenne31::new(0), c1: Mersenne31::new(0) },
        c1: Fp2 { c0: Mersenne31::new(1), c1: Mersenne31::new(0) },
        c2: Fp2 { c0: Mersenne31::new(0), c1: Mersenne31::new(0) },
    };

    // GAMMA^((p-1)/2) in Fp6.
    const FROBENIUS_COEFF: Fp6<Self> = Fp6 {
        c0: Fp2 { c0: Mersenne31::new(953_627_146), c1: Mersenne31::new(879_341_041) },
        c1: Fp2 { c0: Mersenne31::new(0), c1: Mersenne31::new(0) },
        c2: Fp2 { c0: Mersenne31::new(0), c1: Mersenne31::new(0) },
    };
}

pub type Mersenne31Fp2 = Fp2<Mersenne31Config>;
pub type Mersenne31Fp6 = Fp6<Mersenne31Config>;
pub type Mersenne31Fp12 = Fp12<Mersenne31Config>;

#[cfg(test)]
mod ext_tests {
    use super::*;
    use crate::{Mersenne31Backend, FpBackend, Frobenius};
    use proptest::prelude::*;

    fn fm(v: u32) -> Mersenne31 {
        Mersenne31::new(v % Mersenne31Backend::MODULUS)
    }

    fn fe2(c0: u32, c1: u32) -> Mersenne31Fp2 {
        Fp2::new(fm(c0), fm(c1))
    }

    fn canonical2(x: Mersenne31Fp2) -> (u32, u32) {
        (x.c0.value, x.c1.value)
    }

    fn fe6(c0: (u32, u32), c1: (u32, u32), c2: (u32, u32)) -> Mersenne31Fp6 {
        Fp6::new(fe2(c0.0, c0.1), fe2(c1.0, c1.1), fe2(c2.0, c2.1))
    }

    type Coord6 = ((u32, u32), (u32, u32), (u32, u32));

    fn canonical6(x: Mersenne31Fp6) -> Coord6 {
        (canonical2(x.c0), canonical2(x.c1), canonical2(x.c2))
    }

    fn fe12(c0: Coord6, c1: Coord6) -> Mersenne31Fp12 {
        Fp12::new(fe6(c0.0, c0.1, c0.2), fe6(c1.0, c1.1, c1.2))
    }

    fn canonical12(x: Mersenne31Fp12) -> (Coord6, Coord6) {
        (canonical6(x.c0), canonical6(x.c1))
    }

    // Square-and-multiply x^exp using only the ring's own Mul, so it shares
    // no code with the formula-based Frobenius impl and can serve as an
    // independent check on it -- exponents here are p itself (~2^31), so
    // naive repeated multiplication (as field-core's tiny-modulus fixtures
    // use) would be far too slow; this is the same square-and-multiply
    // shape as Fp::pow, just generic over any Copy + Mul ring element.
    fn pow_by_squaring<T: Copy + std::ops::Mul<Output = T>>(base: T, mut exp: u64, one: T) -> T {
        let mut result = one;
        let mut b = base;
        while exp > 0 {
            if exp & 1 == 1 {
                result = result * b;
            }
            b = b * b;
            exp >>= 1;
        }
        result
    }

    const P: u64 = Mersenne31Backend::MODULUS as u64;

    #[test]
    fn frobenius_coeff_matches_euler_criterion() {
        // BETA is a quadratic non-residue mod p, so BETA^((p-1)/2) == -1.
        assert_eq!(
            <Mersenne31Config as Fp2Config>::FROBENIUS_COEFF.value,
            Mersenne31Backend::MODULUS - 1
        );
    }

    #[test]
    fn fp2_frobenius_matches_x_to_the_p() {
        let a = fe2(123_456, 789_012);
        let one = fe2(1, 0);

        assert_eq!(canonical2(a.frobenius()), canonical2(pow_by_squaring(a, P, one)));
    }

    #[test]
    fn fp6_frobenius_matches_x_to_the_p() {
        let a = fe6((1, 2), (3, 4), (5, 6));
        let one = fe6((1, 0), (0, 0), (0, 0));

        assert_eq!(canonical6(a.frobenius()), canonical6(pow_by_squaring(a, P, one)));
    }

    #[test]
    fn fp12_frobenius_matches_x_to_the_p() {
        let a = fe12(((1, 2), (3, 4), (5, 6)), ((7, 8), (9, 10), (11, 12)));
        let one = fe12(((1, 0), (0, 0), (0, 0)), ((0, 0), (0, 0), (0, 0)));

        assert_eq!(canonical12(a.frobenius()), canonical12(pow_by_squaring(a, P, one)));
    }

    proptest! {
        #[test]
        fn fp2_mul_distributes_over_add(a0 in 0u32..1000, a1 in 0u32..1000, b0 in 0u32..1000, b1 in 0u32..1000, c0 in 0u32..1000, c1 in 0u32..1000) {
            let a = fe2(a0, a1);
            let b = fe2(b0, b1);
            let c = fe2(c0, c1);

            let lhs = a * (b + c);
            let rhs = a * b + a * c;

            prop_assert_eq!(canonical2(lhs), canonical2(rhs));
        }

        #[test]
        fn fp2_inverse_times_self_is_one(c0 in 1u32..1000, c1 in 0u32..1000) {
            let a = fe2(c0, c1);
            let product = a * a.inverse();

            prop_assert_eq!(canonical2(product), (1, 0));
        }

        #[test]
        fn fp6_inverse_times_self_is_one(c0 in 1u32..1000, c1 in 0u32..1000, c2 in 0u32..1000) {
            let a = fe6((c0, c1), (c2, 0), (0, 0));
            let product = a * a.inverse();

            prop_assert_eq!(canonical6(product), ((1, 0), (0, 0), (0, 0)));
        }

        #[test]
        fn fp12_inverse_times_self_is_one(a0 in 1u32..1000, a1 in 0u32..1000) {
            let a = fe12(((a0, a1), (0, 0), (0, 0)), ((0, 0), (0, 0), (0, 0)));
            let product = a * a.inverse();

            prop_assert_eq!(canonical12(product), (((1, 0), (0, 0), (0, 0)), ((0, 0), (0, 0), (0, 0))));
        }
    }
}
