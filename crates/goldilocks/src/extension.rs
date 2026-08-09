use field_core::{Fp2, Fp2Config, Fp6, Fp6Config, Fp12, Fp12Config};

use crate::{Goldilocks, GoldilocksConfig};

// Fp2/Fp6/Fp12 tower over Goldilocks. Same "verified offline" convention as
// mersenne31's tower (see its comment for why): GoldilocksBackend stores
// values in plain canonical form, not Montgomery form, so there's no
// to_mont/pow_mont_* compile-time machinery available here -- every
// constant below is the canonical residue found by the accompanying offline
// derivation script (small-non-residue search + modular exponentiation over
// the exact same QuadExt/CubicExt tower arithmetic this crate builds on).
impl Fp2Config for GoldilocksConfig {
    type Base = Goldilocks;

    // 7 is the smallest quadratic non-residue mod p = 2^64-2^32+1.
    const BETA: Self::Base = Goldilocks::new(7);

    // BETA^((p-1)/2) mod p.
    const FROBENIUS_COEFF: Self::Base = Goldilocks::new(18_446_744_069_414_584_320);
}

impl Fp6Config for GoldilocksConfig {
    // u (embedded via Fp2's c1, c0 = 0) is the smallest cubic non-residue in
    // Fp2 found by the search.
    const XI: Fp2<Self> = Fp2 { c0: Goldilocks::new(0), c1: Goldilocks::new(1) };

    // XI^((p-1)/3) and XI^(2(p-1)/3) in Fp2.
    const FROBENIUS_COEFF_C1: Fp2<Self> =
        Fp2 { c0: Goldilocks::new(18_446_744_065_119_617_026), c1: Goldilocks::new(0) };
    const FROBENIUS_COEFF_C2: Fp2<Self> =
        Fp2 { c0: Goldilocks::new(18_446_744_065_119_617_025), c1: Goldilocks::new(0) };
}

impl Fp12Config for GoldilocksConfig {
    // u*v (Fp2's own generator u, embedded via Fp6's c1 with c0 = c2 = 0) is
    // the smallest quadratic non-residue in Fp6 found by the search.
    const GAMMA: Fp6<Self> = Fp6 {
        c0: Fp2 { c0: Goldilocks::new(0), c1: Goldilocks::new(0) },
        c1: Fp2 { c0: Goldilocks::new(0), c1: Goldilocks::new(1) },
        c2: Fp2 { c0: Goldilocks::new(0), c1: Goldilocks::new(0) },
    };

    // GAMMA^((p-1)/2) in Fp6.
    const FROBENIUS_COEFF: Fp6<Self> = Fp6 {
        c0: Fp2 { c0: Goldilocks::new(18_446_744_065_119_617_025), c1: Goldilocks::new(0) },
        c1: Fp2 { c0: Goldilocks::new(0), c1: Goldilocks::new(0) },
        c2: Fp2 { c0: Goldilocks::new(0), c1: Goldilocks::new(0) },
    };
}

pub type GoldilocksFp2 = Fp2<GoldilocksConfig>;
pub type GoldilocksFp6 = Fp6<GoldilocksConfig>;
pub type GoldilocksFp12 = Fp12<GoldilocksConfig>;

#[cfg(test)]
mod ext_tests {
    use super::*;
    use crate::{GoldilocksBackend, FpBackend, Frobenius};
    use proptest::prelude::*;

    fn fg(v: u64) -> Goldilocks {
        Goldilocks::new(v % GoldilocksBackend::MODULUS)
    }

    fn fe2(c0: u64, c1: u64) -> GoldilocksFp2 {
        Fp2::new(fg(c0), fg(c1))
    }

    fn canonical2(x: GoldilocksFp2) -> (u64, u64) {
        (x.c0.value, x.c1.value)
    }

    fn fe6(c0: (u64, u64), c1: (u64, u64), c2: (u64, u64)) -> GoldilocksFp6 {
        Fp6::new(fe2(c0.0, c0.1), fe2(c1.0, c1.1), fe2(c2.0, c2.1))
    }

    type Coord6 = ((u64, u64), (u64, u64), (u64, u64));

    fn canonical6(x: GoldilocksFp6) -> Coord6 {
        (canonical2(x.c0), canonical2(x.c1), canonical2(x.c2))
    }

    fn fe12(c0: Coord6, c1: Coord6) -> GoldilocksFp12 {
        Fp12::new(fe6(c0.0, c0.1, c0.2), fe6(c1.0, c1.1, c1.2))
    }

    fn canonical12(x: GoldilocksFp12) -> (Coord6, Coord6) {
        (canonical6(x.c0), canonical6(x.c1))
    }

    // See mersenne31's ext_tests for why this is square-and-multiply rather
    // than naive repeated multiplication: exponents here are p itself
    // (~2^64), so naive repeated multiplication would be far too slow.
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

    const P: u64 = GoldilocksBackend::MODULUS;

    #[test]
    fn frobenius_coeff_matches_euler_criterion() {
        // BETA is a quadratic non-residue mod p, so BETA^((p-1)/2) == -1.
        assert_eq!(
            <GoldilocksConfig as Fp2Config>::FROBENIUS_COEFF.value,
            GoldilocksBackend::MODULUS - 1
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
        fn fp2_mul_distributes_over_add(a0 in 0u64..1000, a1 in 0u64..1000, b0 in 0u64..1000, b1 in 0u64..1000, c0 in 0u64..1000, c1 in 0u64..1000) {
            let a = fe2(a0, a1);
            let b = fe2(b0, b1);
            let c = fe2(c0, c1);

            let lhs = a * (b + c);
            let rhs = a * b + a * c;

            prop_assert_eq!(canonical2(lhs), canonical2(rhs));
        }

        #[test]
        fn fp2_inverse_times_self_is_one(c0 in 1u64..1000, c1 in 0u64..1000) {
            let a = fe2(c0, c1);
            let product = a * a.inverse();

            prop_assert_eq!(canonical2(product), (1, 0));
        }

        #[test]
        fn fp6_inverse_times_self_is_one(c0 in 1u64..1000, c1 in 0u64..1000, c2 in 0u64..1000) {
            let a = fe6((c0, c1), (c2, 0), (0, 0));
            let product = a * a.inverse();

            prop_assert_eq!(canonical6(product), ((1, 0), (0, 0), (0, 0)));
        }

        #[test]
        fn fp12_inverse_times_self_is_one(a0 in 1u64..1000, a1 in 0u64..1000) {
            let a = fe12(((a0, a1), (0, 0), (0, 0)), ((0, 0), (0, 0), (0, 0)));
            let product = a * a.inverse();

            prop_assert_eq!(canonical12(product), (((1, 0), (0, 0), (0, 0)), ((0, 0), (0, 0), (0, 0))));
        }
    }
}
