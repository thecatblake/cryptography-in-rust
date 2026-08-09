use field_core::{
    Fp, Fp2, Fp2Config, Fp6, Fp6Config, Fp12, Fp12Config, MontFieldConfig, pow_mont_fp2_u32, pow_mont_u32,
    to_mont_u32,
};

use crate::{BabyBear, BabyBearConfig};

// Fp2/Fp6/Fp12 tower over BabyBear. Unlike mersenne31/goldilocks (plain
// canonical backend), BabyBearBackend is Montgomery-form, so BETA and the
// Fp2/Fp6-level Frobenius coefficients can be computed entirely at compile
// time from their canonical (non-Montgomery) values via the to_mont_u32/
// pow_mont_u32/pow_mont_fp2_u32 helpers -- the same pattern field-core's own
// Mod17Mont/Mod13Mont test fixtures use. GAMMA and Fp12's FROBENIUS_COEFF
// are the one exception: no pow_mont_fp6 helper exists (see Fp12Config's
// doc comment in field-core), so those two are still the canonical values
// found by the accompanying offline derivation script, converted to
// Montgomery form the same to_mont_u32 way.
impl Fp2Config for BabyBearConfig {
    type Base = BabyBear;

    // 11 is the smallest quadratic non-residue mod p.
    const BETA: Self::Base =
        Fp::new(to_mont_u32(11, Self::R2, Self::N_PRIME, Self::MODULUS));

    const FROBENIUS_COEFF: Self::Base = Fp::new(pow_mont_u32(
        <Self as Fp2Config>::BETA.value,
        (Self::MODULUS - 1) / 2,
        Self::R2,
        Self::N_PRIME,
        Self::MODULUS,
    ));
}

impl Fp6Config for BabyBearConfig {
    // 1 + u is the smallest cubic non-residue in Fp2 found by the search.
    const XI: Fp2<Self> = Fp2 {
        c0: Fp::new(to_mont_u32(1, Self::R2, Self::N_PRIME, Self::MODULUS)),
        c1: Fp::new(to_mont_u32(1, Self::R2, Self::N_PRIME, Self::MODULUS)),
    };

    const FROBENIUS_COEFF_C1: Fp2<Self> = {
        let (c0, c1) = pow_mont_fp2_u32(
            <Self as Fp6Config>::XI.c0.value,
            <Self as Fp6Config>::XI.c1.value,
            (Self::MODULUS - 1) / 3,
            <Self as Fp2Config>::BETA.value,
            Self::R2,
            Self::N_PRIME,
            Self::MODULUS,
        );
        Fp2 { c0: Fp::new(c0), c1: Fp::new(c1) }
    };
    const FROBENIUS_COEFF_C2: Fp2<Self> = {
        let (c0, c1) = pow_mont_fp2_u32(
            <Self as Fp6Config>::XI.c0.value,
            <Self as Fp6Config>::XI.c1.value,
            2 * ((Self::MODULUS - 1) / 3),
            <Self as Fp2Config>::BETA.value,
            Self::R2,
            Self::N_PRIME,
            Self::MODULUS,
        );
        Fp2 { c0: Fp::new(c0), c1: Fp::new(c1) }
    };
}

// Canonical (non-Montgomery) coordinates of BabyBearConfig::GAMMA -- v
// (Fp6's own generator, embedded via c1 with c0 = c2 = 0), the smallest
// quadratic non-residue in Fp6 found by the offline derivation script.
const GAMMA_CANONICAL: ((u32, u32), (u32, u32), (u32, u32)) = ((0, 0), (1, 0), (0, 0));

// Canonical coordinates of BabyBearConfig::FROBENIUS_COEFF = GAMMA^((p-1)/2)
// in Fp6, computed by the same offline derivation script and cross-checked
// there against x^p via repeated Fp12 multiplication (see the script for
// details -- the same "verified offline" convention field-core's own Fp12
// fixtures use, since no pow_mont_fp6 helper exists).
const FROBENIUS_COEFF_CANONICAL: ((u32, u32), (u32, u32), (u32, u32)) =
    ((1_580_496_498, 1_982_416_194), (0, 0), (0, 0));

impl Fp12Config for BabyBearConfig {
    const GAMMA: Fp6<Self> = Fp6 {
        c0: Fp2 {
            c0: Fp::new(to_mont_u32(GAMMA_CANONICAL.0.0, Self::R2, Self::N_PRIME, Self::MODULUS)),
            c1: Fp::new(to_mont_u32(GAMMA_CANONICAL.0.1, Self::R2, Self::N_PRIME, Self::MODULUS)),
        },
        c1: Fp2 {
            c0: Fp::new(to_mont_u32(GAMMA_CANONICAL.1.0, Self::R2, Self::N_PRIME, Self::MODULUS)),
            c1: Fp::new(to_mont_u32(GAMMA_CANONICAL.1.1, Self::R2, Self::N_PRIME, Self::MODULUS)),
        },
        c2: Fp2 {
            c0: Fp::new(to_mont_u32(GAMMA_CANONICAL.2.0, Self::R2, Self::N_PRIME, Self::MODULUS)),
            c1: Fp::new(to_mont_u32(GAMMA_CANONICAL.2.1, Self::R2, Self::N_PRIME, Self::MODULUS)),
        },
    };

    const FROBENIUS_COEFF: Fp6<Self> = Fp6 {
        c0: Fp2 {
            c0: Fp::new(to_mont_u32(
                FROBENIUS_COEFF_CANONICAL.0.0,
                Self::R2,
                Self::N_PRIME,
                Self::MODULUS,
            )),
            c1: Fp::new(to_mont_u32(
                FROBENIUS_COEFF_CANONICAL.0.1,
                Self::R2,
                Self::N_PRIME,
                Self::MODULUS,
            )),
        },
        c1: Fp2 {
            c0: Fp::new(to_mont_u32(
                FROBENIUS_COEFF_CANONICAL.1.0,
                Self::R2,
                Self::N_PRIME,
                Self::MODULUS,
            )),
            c1: Fp::new(to_mont_u32(
                FROBENIUS_COEFF_CANONICAL.1.1,
                Self::R2,
                Self::N_PRIME,
                Self::MODULUS,
            )),
        },
        c2: Fp2 {
            c0: Fp::new(to_mont_u32(
                FROBENIUS_COEFF_CANONICAL.2.0,
                Self::R2,
                Self::N_PRIME,
                Self::MODULUS,
            )),
            c1: Fp::new(to_mont_u32(
                FROBENIUS_COEFF_CANONICAL.2.1,
                Self::R2,
                Self::N_PRIME,
                Self::MODULUS,
            )),
        },
    };
}

pub type BabyBearFp2 = Fp2<BabyBearConfig>;
pub type BabyBearFp6 = Fp6<BabyBearConfig>;
pub type BabyBearFp12 = Fp12<BabyBearConfig>;

#[cfg(test)]
mod ext_tests {
    use super::*;
    use crate::{BabyBearBackend, FpBackend, Frobenius};
    use proptest::prelude::*;

    fn fb(v: u32) -> BabyBear {
        BabyBear::new(BabyBearBackend::to_mont(v % BabyBearBackend::MODULUS))
    }

    fn canonical(x: BabyBear) -> u32 {
        BabyBearBackend::from_mont(x.value)
    }

    fn fe2(c0: u32, c1: u32) -> BabyBearFp2 {
        Fp2::new(fb(c0), fb(c1))
    }

    fn canonical2(x: BabyBearFp2) -> (u32, u32) {
        (canonical(x.c0), canonical(x.c1))
    }

    fn fe6(c0: (u32, u32), c1: (u32, u32), c2: (u32, u32)) -> BabyBearFp6 {
        Fp6::new(fe2(c0.0, c0.1), fe2(c1.0, c1.1), fe2(c2.0, c2.1))
    }

    type Coord6 = ((u32, u32), (u32, u32), (u32, u32));

    fn canonical6(x: BabyBearFp6) -> Coord6 {
        (canonical2(x.c0), canonical2(x.c1), canonical2(x.c2))
    }

    fn fe12(c0: Coord6, c1: Coord6) -> BabyBearFp12 {
        Fp12::new(fe6(c0.0, c0.1, c0.2), fe6(c1.0, c1.1, c1.2))
    }

    fn canonical12(x: BabyBearFp12) -> (Coord6, Coord6) {
        (canonical6(x.c0), canonical6(x.c1))
    }

    // See mersenne31's ext_tests for why this is square-and-multiply rather
    // than naive repeated multiplication. Works directly on Montgomery-form
    // values: Montgomery multiplication is a ring isomorphism, so repeated
    // squaring on Montgomery-form operands computes (x^exp)_mont, as long as
    // `one` is Montgomery-form 1 too (matching Fp<MontWideBackend<_>>::pow).
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

    const P: u64 = BabyBearBackend::MODULUS as u64;

    #[test]
    fn frobenius_coeff_matches_euler_criterion() {
        // BETA is a quadratic non-residue mod p, so BETA^((p-1)/2) == -1.
        assert_eq!(canonical(<BabyBearConfig as Fp2Config>::FROBENIUS_COEFF), BabyBearBackend::MODULUS - 1);
    }

    #[test]
    fn fp2_frobenius_matches_x_to_the_p() {
        let a = fe2(123_456, 78_901);
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
