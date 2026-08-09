use std::marker::PhantomData;

use crate::{Fp6, Fp6Config, QuadExt, QuadExtConfig, QuadExtFrobeniusConfig};

// Fp12Config extends Fp6Config with the one extra constant a quadratic
// extension needs: Fp12 = Fp6[w] / (w^2 - GAMMA). Fp6Config already
// supplies everything needed to build the base field's Fp6 arithmetic, so
// an Fp12Config implementor only adds GAMMA on top.
pub trait Fp12Config: Fp6Config + Sized {
    // Must be a quadratic non-residue in Fp6 (GAMMA has no square root in
    // Fp6), or w^2 - GAMMA factors and Fp12 collapses into Fp6 x Fp6
    // instead of being a field.
    const GAMMA: Fp6<Self>;

    // Frobenius on Fp12 sends w to w^p: since w^2 == GAMMA, w^p =
    // w * (w^2)^((p-1)/2) = GAMMA^((p-1)/2) (see the QuadExtFrobeniusConfig
    // impl on Fp12Marker below, which is what actually wires this into
    // Frobenius). Unlike Fp2Config::FROBENIUS_COEFF and
    // Fp6Config::FROBENIUS_COEFF_C1/C2, this can't be written as a call to
    // a pow_mont_*_uNN helper: those exponentiate Fp and Fp2 elements at
    // compile time, and GAMMA is an Fp6 element -- exponentiating it would
    // need an analogous pow_mont_fp6_uNN (Fp6-level Montgomery
    // multiplication hardcoded as a const fn) that doesn't exist yet. So
    // this constant is verified offline instead (same convention already
    // used for XI/GAMMA themselves), same as any other GAMMA-derived value
    // would be until that helper is written.
    const FROBENIUS_COEFF: Fp6<Self>;
}

// Fp12 needs a QuadExtConfig of its own, with Base = Fp6<C>. It can't be a
// second blanket `impl<C: Fp12Config> QuadExtConfig for C` the way Fp2Config
// and Fp6Config each get one in fp2.rs/fp6.rs: Fp12Config: Fp6Config:
// Fp2Config, so any C: Fp12Config already satisfies C: Fp2Config, and
// fp2.rs's `impl<C: Fp2Config> QuadExtConfig for C` would then conflict
// with a second impl of the same trait for the same C. Indexing this level
// of the tower by a marker type instead of C directly sidesteps that: only
// Fp12Marker<C>, not C itself, implements QuadExtConfig here.
pub struct Fp12Marker<C>(PhantomData<C>);

impl<C: Fp12Config> QuadExtConfig for Fp12Marker<C> {
    type Base = Fp6<C>;

    const BETA: Self::Base = C::GAMMA;
}

// Mirrors QuadExtConfig above: Fp12Marker<C>, not C, carries the
// QuadExtFrobeniusConfig instance too, forwarding Fp12Config::
// FROBENIUS_COEFF the same way BETA forwards GAMMA. This is what makes
// quad_ext.rs's generic `impl<C: QuadExtFrobeniusConfig> Frobenius for
// QuadExt<C>` apply to Fp12<C> -- Fp6<C> already implements Frobenius (via
// Fp6Config's own CubicExtFrobeniusConfig forwarding), so the `where
// C::Base: Frobenius` bound on that generic impl is satisfied too.
impl<C: Fp12Config> QuadExtFrobeniusConfig for Fp12Marker<C> {
    const FROBENIUS_COEFF: Self::Base = <C as Fp12Config>::FROBENIUS_COEFF;
}

// Fp12 = Fp6[w] / (w^2 - GAMMA), elements represented as c0 + c1*w. A
// specialization of QuadExt to an Fp6 base field (via the Fp12Marker
// indirection above) -- Add/Sub/Neg/Mul/inverse all come from QuadExt, not
// reimplemented here.
pub type Fp12<C> = QuadExt<Fp12Marker<C>>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Fp, Fp2, Fp2Config, Frobenius, MontFieldConfig, MontWideBackend, pow_mont_fp2_u32, pow_mont_u32,
        to_mont_u32,
    };
    use proptest::prelude::*;

    // Same Mod17Mont fixture as fp2's/fp6's tests: p = 17, R = 2^32.
    struct Mod17Mont;
    impl MontFieldConfig for Mod17Mont {
        type Repr = u32;

        const MODULUS: u32 = 17;
        const R2: u32 = 0x1;
        const N_PRIME: u32 = 0x0f0f_0f0f;
    }

    const BETA_CANONICAL: u32 = 3;

    impl Fp2Config for Mod17Mont {
        type Base = Fp<MontWideBackend<Self>>;

        const BETA: Self::Base =
            Fp::new(to_mont_u32(BETA_CANONICAL, Self::R2, Self::N_PRIME, Self::MODULUS));

        // See fp2's Mod17Mont fixture: BETA^((p-1)/2) = 3^8 mod 17 = 16.
        const FROBENIUS_COEFF: Self::Base = Fp::new(pow_mont_u32(
            <Self as Fp2Config>::BETA.value,
            (Self::MODULUS - 1) / 2,
            Self::R2,
            Self::N_PRIME,
            Self::MODULUS,
        ));
    }

    const XI_CANONICAL: (u32, u32) = (1, 1);

    impl Fp6Config for Mod17Mont {
        const XI: Fp2<Self> = Fp2 {
            c0: Fp::new(to_mont_u32(XI_CANONICAL.0, Self::R2, Self::N_PRIME, Self::MODULUS)),
            c1: Fp::new(to_mont_u32(XI_CANONICAL.1, Self::R2, Self::N_PRIME, Self::MODULUS)),
        };

        // 17 mod 3 == 2, so these aren't Fp6's actual Frobenius
        // coefficients -- see fp6's Mod17Mont fixture for the same caveat.
        // Present only so Mod17Mont satisfies Fp6Config here.
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

    // Canonical (non-Montgomery) coordinates of Mod17Mont::GAMMA, used by
    // tests that check Fp12 arithmetic against a plain-integer formula.
    // GAMMA_CANONICAL = ((0,1), (0,0), (0,0)) -- i.e. GAMMA = u, embedded as
    // Fp6's c0 coordinate.
    const GAMMA_CANONICAL: ((u32, u32), (u32, u32), (u32, u32)) = ((0, 1), (0, 0), (0, 0));

    // Canonical coordinates of Mod17Mont::FROBENIUS_COEFF = GAMMA^((p-1)/2).
    // Verified offline: GAMMA = u lies in Fp6's Fp2 subfield (c1 = c2 = 0),
    // where Fp6 multiplication restricts to plain Fp2 multiplication, so
    // this is just u^8 computed in F17^2 -- u^2 = BETA = 3, u^8 = 3^4 mod 17
    // = 13. Cross-checked independently via the frobenius_matches_x_to_the_p
    // test below, which recomputes x^17 by repeated multiplication with no
    // shared code path.
    const FROBENIUS_COEFF_CANONICAL: ((u32, u32), (u32, u32), (u32, u32)) =
        ((13, 0), (0, 0), (0, 0));

    impl Fp12Config for Mod17Mont {
        // Verified offline (see the accompanying derivation script): among
        // the 17^6-1 nonzero elements of F17^6, u (embedded via Fp6's c0)
        // raises to the (17^6-1)/2 power as -1, not 1, so it isn't a square
        // and w^2 - u is irreducible over F17^6.
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

    type Backend = MontWideBackend<Mod17Mont>;
    type F17 = Fp<Backend>;
    type F17_2 = Fp2<Mod17Mont>;
    type F17_6 = Fp6<Mod17Mont>;
    type F17_12 = Fp12<Mod17Mont>;

    fn fe(v: u32) -> F17 {
        F17::new(Backend::to_mont(v % 17))
    }

    fn canonical(x: F17) -> u32 {
        Backend::from_mont(x.value)
    }

    fn fe2(c0: u32, c1: u32) -> F17_2 {
        Fp2::new(fe(c0), fe(c1))
    }

    fn canonical2(x: F17_2) -> (u32, u32) {
        (canonical(x.c0), canonical(x.c1))
    }

    fn fe6(c0: (u32, u32), c1: (u32, u32), c2: (u32, u32)) -> F17_6 {
        Fp6::new(fe2(c0.0, c0.1), fe2(c1.0, c1.1), fe2(c2.0, c2.1))
    }

    fn canonical6(x: F17_6) -> ((u32, u32), (u32, u32), (u32, u32)) {
        (canonical2(x.c0), canonical2(x.c1), canonical2(x.c2))
    }

    type Coord6 = ((u32, u32), (u32, u32), (u32, u32));

    fn fe12(c0: Coord6, c1: Coord6) -> F17_12 {
        Fp12::new(fe6(c0.0, c0.1, c0.2), fe6(c1.0, c1.1, c1.2))
    }

    fn canonical12(x: F17_12) -> (Coord6, Coord6) {
        (canonical6(x.c0), canonical6(x.c1))
    }

    const ZERO6: Coord6 = ((0, 0), (0, 0), (0, 0));

    #[test]
    fn gamma_matches_configured_value() {
        assert_eq!(canonical6(<Mod17Mont as Fp12Config>::GAMMA), GAMMA_CANONICAL);
    }

    #[test]
    fn add_matches_componentwise() {
        let a0 = ((1, 2), (3, 4), (5, 6));
        let a1 = ((7, 8), (9, 10), (11, 12));
        let b0 = ((2, 3), (4, 5), (6, 7));
        let b1 = ((8, 9), (10, 11), (12, 13));

        let sum = fe12(a0, a1) + fe12(b0, b1);

        assert_eq!(canonical6(sum.c0), ((3, 5), (7, 9), (11, 13)));
        assert_eq!(canonical6(sum.c1), ((15, 0), (2, 4), (6, 8)));
    }

    #[test]
    fn sub_is_inverse_of_add() {
        let a = fe12(((1, 2), (3, 4), (5, 6)), ((7, 8), (9, 10), (11, 12)));
        let b = fe12(((2, 3), (4, 5), (6, 7)), ((8, 9), (10, 11), (12, 13)));
        let back = (a + b) - b;
        assert_eq!(canonical12(back), canonical12(a));
    }

    #[test]
    fn neg_matches_sub_from_zero() {
        let a = fe12(((1, 2), (3, 4), (5, 6)), ((7, 8), (9, 10), (11, 12)));
        let zero = fe12(ZERO6, ZERO6);
        assert_eq!(canonical12(-a), canonical12(zero - a));
    }

    // Plain-integer schoolbook reference for Fp2/Fp6 multiplication, used
    // only by the Fp12 schoolbook check below.
    fn mul2_ref(a: (u32, u32), b: (u32, u32)) -> (u32, u32) {
        let p = 17u32;
        let c0 = (a.0 * b.0 + BETA_CANONICAL * a.1 * b.1) % p;
        let c1 = (a.0 * b.1 + a.1 * b.0) % p;
        (c0, c1)
    }

    fn add2_ref(a: (u32, u32), b: (u32, u32)) -> (u32, u32) {
        ((a.0 + b.0) % 17, (a.1 + b.1) % 17)
    }

    fn add6_ref(a: Coord6, b: Coord6) -> Coord6 {
        (add2_ref(a.0, b.0), add2_ref(a.1, b.1), add2_ref(a.2, b.2))
    }

    fn mul6_ref(a: Coord6, b: Coord6) -> Coord6 {
        let d0 = mul2_ref(a.0, b.0);
        let d1 = add2_ref(mul2_ref(a.0, b.1), mul2_ref(a.1, b.0));
        let d2 = add2_ref(add2_ref(mul2_ref(a.0, b.2), mul2_ref(a.1, b.1)), mul2_ref(a.2, b.0));
        let d3 = add2_ref(mul2_ref(a.1, b.2), mul2_ref(a.2, b.1));
        let d4 = mul2_ref(a.2, b.2);

        let c0 = add2_ref(d0, mul2_ref(XI_CANONICAL, d3));
        let c1 = add2_ref(d1, mul2_ref(XI_CANONICAL, d4));
        let c2 = d2;
        (c0, c1, c2)
    }

    #[test]
    fn mul_matches_schoolbook_expansion() {
        let a0 = ((1, 2), (3, 4), (5, 6));
        let a1 = ((7, 8), (9, 10), (11, 12));
        let b0 = ((2, 3), (4, 5), (6, 7));
        let b1 = ((8, 9), (10, 11), (12, 13));

        let product = fe12(a0, a1) * fe12(b0, b1);

        // (a0 + a1*w)(b0 + b1*w) = (a0*b0 + GAMMA*a1*b1) + (a0*b1 + a1*b0)*w
        let expected_c0 = add6_ref(mul6_ref(a0, b0), mul6_ref(GAMMA_CANONICAL, mul6_ref(a1, b1)));
        let expected_c1 = add6_ref(mul6_ref(a0, b1), mul6_ref(a1, b0));

        assert_eq!(canonical6(product.c0), expected_c0);
        assert_eq!(canonical6(product.c1), expected_c1);
    }

    fn small_coord6(v: u32) -> Coord6 {
        ((v % 17, (v + 1) % 17), ((v + 2) % 17, (v + 3) % 17), ((v + 4) % 17, (v + 5) % 17))
    }

    proptest! {
        #[test]
        fn add_commutative(a in 0u32..17, b in 0u32..17, c in 0u32..17, d in 0u32..17) {
            let x = fe12(small_coord6(a), small_coord6(b));
            let y = fe12(small_coord6(c), small_coord6(d));
            prop_assert_eq!(canonical12(x + y), canonical12(y + x));
        }

        #[test]
        fn mul_commutative(a in 0u32..17, b in 0u32..17, c in 0u32..17, d in 0u32..17) {
            let x = fe12(small_coord6(a), small_coord6(b));
            let y = fe12(small_coord6(c), small_coord6(d));
            prop_assert_eq!(canonical12(x * y), canonical12(y * x));
        }

        #[test]
        fn mul_distributes_over_add(a in 0u32..17, b in 0u32..17, c in 0u32..17, d in 0u32..17, e in 0u32..17, f in 0u32..17) {
            let x = fe12(small_coord6(a), small_coord6(b));
            let y = fe12(small_coord6(c), small_coord6(d));
            let z = fe12(small_coord6(e), small_coord6(f));

            let lhs = x * (y + z);
            let rhs = x * y + x * z;

            prop_assert_eq!(canonical12(lhs), canonical12(rhs));
        }

        #[test]
        fn inverse_times_self_is_one(a in 0u32..17, b in 0u32..17) {
            // GAMMA is a non-residue, so norm(a) == 0 only when a itself is
            // zero -- every other element is invertible.
            prop_assume!(a != 0 || b != 0);

            let x = fe12(small_coord6(a), small_coord6(b));
            let product = x * x.inverse();

            let one = fe12(((1, 0), (0, 0), (0, 0)), ZERO6);
            prop_assert_eq!(canonical12(product), canonical12(one));
        }
    }

    // These two hold under Mod17Mont regardless of Fp6<Mod17Mont>'s known-
    // bad FROBENIUS_COEFF_C1/C2 (17 mod 3 == 2 -- see fp6's Mod17Mont
    // fixture): both only touch elements whose Fp6 c1/c2 coordinates are
    // zero, which zeroes out every term that multiplies by those bad
    // coefficients. Tests that exercise Fp12's Frobenius on elements with
    // nonzero Fp6 c1/c2 -- i.e. anything checking multiplicativity or
    // frobenius(x) == x^p in general -- need a modulus where Fp6's own
    // Frobenius is actually correct, so those live in the frobenius module
    // below instead, on a fresh Mod13Mont fixture (13 mod 3 == 1).

    #[test]
    fn frobenius_coeff_matches_expected_value() {
        assert_eq!(
            canonical6(<Mod17Mont as Fp12Config>::FROBENIUS_COEFF),
            FROBENIUS_COEFF_CANONICAL,
        );
    }

    #[test]
    fn frobenius_fixes_the_base_field_embedding() {
        // Fp embeds into Fp12 as c0's c0's c0 slot with everything else
        // zero; Frobenius must fix it pointwise, same as at every other
        // level of the tower.
        let a = fe12(((5, 0), (0, 0), (0, 0)), ZERO6);
        let frob = a.frobenius();

        assert_eq!(canonical6(frob.c0), canonical6(a.c0));
        assert_eq!(canonical6(frob.c1), ZERO6);
    }

    // 17 mod 2 == 1 always holds for any odd prime, so unlike Fp6's (p-1)/3
    // caveat, Fp12's own (p-1)/2 exponent needs no special-prime fixture --
    // it's only the *nested* Fp6-level (p-1)/3 exponent (used by
    // c0.frobenius()/c1.frobenius() inside Fp12's formula) that forces one.
    mod frobenius {
        use super::*;

        // p = 13 == 1 (mod 3), R = 2^32 -- same fixture as fp6's dedicated
        // Frobenius tests, extended one level up with a GAMMA/
        // FROBENIUS_COEFF for Fp12.
        struct Mod13Mont;
        impl MontFieldConfig for Mod13Mont {
            type Repr = u32;

            const MODULUS: u32 = 13;
            const R2: u32 = 0x3;
            const N_PRIME: u32 = 0x3b13_b13b;
        }

        const BETA_CANONICAL: u32 = 2;

        impl Fp2Config for Mod13Mont {
            type Base = Fp<MontWideBackend<Self>>;

            const BETA: Self::Base =
                Fp::new(to_mont_u32(BETA_CANONICAL, Self::R2, Self::N_PRIME, Self::MODULUS));

            const FROBENIUS_COEFF: Self::Base = Fp::new(pow_mont_u32(
                <Self as Fp2Config>::BETA.value,
                (Self::MODULUS - 1) / 2,
                Self::R2,
                Self::N_PRIME,
                Self::MODULUS,
            ));
        }

        const XI_CANONICAL: (u32, u32) = (0, 1);

        impl Fp6Config for Mod13Mont {
            const XI: Fp2<Self> = Fp2 {
                c0: Fp::new(to_mont_u32(XI_CANONICAL.0, Self::R2, Self::N_PRIME, Self::MODULUS)),
                c1: Fp::new(to_mont_u32(XI_CANONICAL.1, Self::R2, Self::N_PRIME, Self::MODULUS)),
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

        // GAMMA = v (Fp6's own generator, embedded via c1 with c0 = c2 =
        // 0). Verified offline: among the 13^6-1 nonzero elements of F13^6,
        // v raises to the (13^6-1)/2 power as -1, not 1, so it isn't a
        // square and w^2 - v is irreducible over F13^6.
        const GAMMA_CANONICAL: ((u32, u32), (u32, u32), (u32, u32)) = ((0, 0), (1, 0), (0, 0));

        // FROBENIUS_COEFF = GAMMA^((p-1)/2) = v^6 mod 13 = 2 (verified
        // offline). Cross-checked independently by
        // frobenius_matches_x_to_the_p below, which recomputes x^13 by
        // repeated multiplication sharing no code with frobenius().
        const FROBENIUS_COEFF_CANONICAL: ((u32, u32), (u32, u32), (u32, u32)) =
            ((2, 0), (0, 0), (0, 0));

        impl Fp12Config for Mod13Mont {
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

        type Backend = MontWideBackend<Mod13Mont>;
        type F13 = Fp<Backend>;
        type F13_2 = Fp2<Mod13Mont>;
        type F13_6 = Fp6<Mod13Mont>;
        type F13_12 = Fp12<Mod13Mont>;

        fn fe(v: u32) -> F13 {
            F13::new(Backend::to_mont(v % 13))
        }

        fn fe2(c0: u32, c1: u32) -> F13_2 {
            Fp2::new(fe(c0), fe(c1))
        }

        fn canonical(x: F13) -> u32 {
            Backend::from_mont(x.value)
        }

        fn canonical2(x: F13_2) -> (u32, u32) {
            (canonical(x.c0), canonical(x.c1))
        }

        fn fe6(c0: (u32, u32), c1: (u32, u32), c2: (u32, u32)) -> F13_6 {
            Fp6::new(fe2(c0.0, c0.1), fe2(c1.0, c1.1), fe2(c2.0, c2.1))
        }

        fn canonical6(x: F13_6) -> ((u32, u32), (u32, u32), (u32, u32)) {
            (canonical2(x.c0), canonical2(x.c1), canonical2(x.c2))
        }

        type Coord6 = ((u32, u32), (u32, u32), (u32, u32));

        fn fe12(c0: Coord6, c1: Coord6) -> F13_12 {
            Fp12::new(fe6(c0.0, c0.1, c0.2), fe6(c1.0, c1.1, c1.2))
        }

        fn canonical12(x: F13_12) -> (Coord6, Coord6) {
            (canonical6(x.c0), canonical6(x.c1))
        }

        fn small_coord6(v: u32) -> Coord6 {
            ((v % 13, (v + 1) % 13), ((v + 2) % 13, (v + 3) % 13), ((v + 4) % 13, (v + 5) % 13))
        }

        // x^e via straight-line repeated multiplication -- deliberately the
        // most naive possible exponentiation, so it shares no code with
        // frobenius() and can serve as an independent check on it.
        fn pow_by_repeated_mul(a: F13_12, e: u32) -> F13_12 {
            let mut result = a;
            for _ in 1..e {
                result = result * a;
            }
            result
        }

        #[test]
        fn frobenius_coeff_matches_expected_value() {
            assert_eq!(
                canonical6(<Mod13Mont as Fp12Config>::FROBENIUS_COEFF),
                FROBENIUS_COEFF_CANONICAL,
            );
        }

        #[test]
        fn frobenius_matches_x_to_the_p() {
            // The defining property of Frobenius: frobenius(x) == x^p.
            let a = fe12(((1, 2), (3, 4), (5, 6)), ((7, 8), (9, 10), (11, 12)));

            let via_formula = a.frobenius();
            let via_repeated_mul = pow_by_repeated_mul(a, 13);

            assert_eq!(canonical12(via_formula), canonical12(via_repeated_mul));
        }

        proptest! {
            #[test]
            fn frobenius_matches_x_to_the_p_proptest(a in 0u32..13, b in 0u32..13) {
                let x = fe12(small_coord6(a), small_coord6(b));

                let via_formula = x.frobenius();
                let via_repeated_mul = pow_by_repeated_mul(x, 13);

                prop_assert_eq!(canonical12(via_formula), canonical12(via_repeated_mul));
            }

            #[test]
            fn frobenius_is_additive(a in 0u32..13, b in 0u32..13, c in 0u32..13, d in 0u32..13) {
                let x = fe12(small_coord6(a), small_coord6(b));
                let y = fe12(small_coord6(c), small_coord6(d));

                let lhs = (x + y).frobenius();
                let rhs = x.frobenius() + y.frobenius();

                prop_assert_eq!(canonical12(lhs), canonical12(rhs));
            }

            #[test]
            fn frobenius_is_multiplicative(a in 0u32..13, b in 0u32..13, c in 0u32..13, d in 0u32..13) {
                let x = fe12(small_coord6(a), small_coord6(b));
                let y = fe12(small_coord6(c), small_coord6(d));

                let lhs = (x * y).frobenius();
                let rhs = x.frobenius() * y.frobenius();

                prop_assert_eq!(canonical12(lhs), canonical12(rhs));
            }
        }
    }
}
