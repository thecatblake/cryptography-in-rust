use crate::{Fp, Frobenius, MontFieldConfig, MontWideBackend, QuadExt, QuadExtConfig};

// Fp2Config extends MontFieldConfig with the one extra constant a quadratic
// extension needs: Fp2 = Fp[u] / (u^2 - BETA). The MontFieldConfig
// supertrait already supplies everything needed to build the base field's
// backend (MontWideBackend<Self>, via MODULUS/R2/N_PRIME), so an Fp2Config
// implementor holds both the base field's Montgomery backend and BETA.
//
// BETA is typed as an Fp element (not a raw Self::Repr) and fixed at
// compile time via Fp::new being a const fn, so it's a genuine associated
// constant -- no conversion or wrapping happens at the use site, same as
// referencing MontFieldConfig::MODULUS/R2/N_PRIME directly.
pub trait Fp2Config: MontFieldConfig + Sized {
    // Must be a quadratic non-residue mod Self::MODULUS, or u^2 - BETA
    // factors and Fp2 collapses to Fp x Fp instead of being a field.
    //
    // Write this as Fp::new(to_mont_uNN(plain_value, Self::R2, Self::N_PRIME,
    // Self::MODULUS)) (see field_core::to_mont_u32 et al.): the plain,
    // canonical residue goes in, and the Montgomery-form conversion runs at
    // compile time as part of evaluating this const, not at runtime.
    const BETA: Fp<MontWideBackend<Self>>;

    // Frobenius on Fp2 sends u to u^p: since u^2 == BETA, u^p =
    // u * (u^2)^((p-1)/2) = FROBENIUS_COEFF * u, so
    // frobenius(c0 + c1*u) = c0 + (FROBENIUS_COEFF*c1)*u (see the Frobenius
    // impl below).
    //
    // Write this as Fp::new(pow_mont_uNN(<Self as Fp2Config>::BETA.value,
    // (Self::MODULUS - 1) / 2, Self::R2, Self::N_PRIME, Self::MODULUS))
    // (see field_core::pow_mont_u32 et al.): BETA is already in Montgomery
    // form, and pow_mont_uNN's square-and-multiply keeps the result there
    // too, so the whole exponentiation runs at compile time as part of
    // evaluating this const, not at runtime.
    const FROBENIUS_COEFF: Fp<MontWideBackend<Self>>;
}

// Every Fp2Config is a QuadExtConfig with Base fixed to the Montgomery
// backend's Fp -- this is what lets Fp2<C> below be a plain alias for
// QuadExt<C> instead of a hand-rolled struct with its own arithmetic impls.
impl<C: Fp2Config> QuadExtConfig for C {
    type Base = Fp<MontWideBackend<C>>;

    const BETA: Self::Base = <C as Fp2Config>::BETA;
}

// Fp2 = Fp[u] / (u^2 - BETA), elements represented as c0 + c1*u. A
// specialization of QuadExt to a Montgomery-backend base field, via
// Fp2Config's blanket QuadExtConfig impl above -- Add/Sub/Neg/Mul/inverse
// all come from QuadExt, not reimplemented here.
pub type Fp2<C> = QuadExt<C>;

// Frobenius on Fp2 = Fp[u]/(u^2-BETA): x^p for x = c0 + c1*u expands, using
// that Frobenius is additive and fixes Fp itself (Fp's Frobenius impl is
// the identity, by Fermat's little theorem), to c0^p + c1^p*u^p = c0 +
// c1*u^p. u^p reduces to FROBENIUS_COEFF*u by Fp2Config::FROBENIUS_COEFF's
// definition, so this is c0 + (FROBENIUS_COEFF*c1)*u -- one base-field
// multiply, no exponentiation at runtime.
impl<C: Fp2Config> Frobenius for Fp2<C> {
    fn frobenius(self) -> Self {
        Fp2::new(self.c0, C::FROBENIUS_COEFF * self.c1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{pow_mont_u32, to_mont_u32};
    use proptest::prelude::*;

    // Same shape as field's Mod17Mont test fixture, just u32-native instead
    // of U256: p = 17, R = 2^32.
    struct Mod17Mont;
    impl MontFieldConfig for Mod17Mont {
        type Repr = u32;

        const MODULUS: u32 = 17;
        // R2 = R^2 mod p, N_PRIME = -p^-1 mod R, computed offline (see
        // field's Mod17Mont for the derivation approach) and round-trip
        // verified by mont_roundtrip_is_identity below.
        const R2: u32 = 0x1;
        const N_PRIME: u32 = 0x0f0f_0f0f;
    }

    // Canonical (non-Montgomery) residue of Mod17Mont::BETA, used by tests
    // that need to check Fp2 arithmetic against a plain-integer formula.
    const BETA_CANONICAL: u32 = 3;

    impl Fp2Config for Mod17Mont {
        // QRs mod 17 are {1,2,4,8,9,13,15,16}; 3 isn't among them, so
        // u^2 - 3 is irreducible over F17. BETA_CANONICAL (the plain
        // residue) goes in; to_mont_u32 converts it to Montgomery form as
        // part of evaluating this const, at compile time.
        const BETA: Fp<MontWideBackend<Self>> =
            Fp::new(to_mont_u32(BETA_CANONICAL, Self::R2, Self::N_PRIME, Self::MODULUS));

        // t = BETA^((p-1)/2) = 3^8 mod 17 = 16 = -1 mod 17 -- expected by
        // Euler's criterion, since BETA is a quadratic non-residue.
        const FROBENIUS_COEFF: Fp<MontWideBackend<Self>> = Fp::new(pow_mont_u32(
            <Self as Fp2Config>::BETA.value,
            (Self::MODULUS - 1) / 2,
            Self::R2,
            Self::N_PRIME,
            Self::MODULUS,
        ));
    }

    type Backend = MontWideBackend<Mod17Mont>;
    type F17 = Fp<Backend>;
    type F17_2 = Fp2<Mod17Mont>;

    fn fe(v: u32) -> F17 {
        F17::new(Backend::to_mont(v % 17))
    }

    fn canonical(x: F17) -> u32 {
        Backend::from_mont(x.value)
    }

    fn fe2(c0: u32, c1: u32) -> F17_2 {
        Fp2::new(fe(c0), fe(c1))
    }

    #[test]
    fn mont_roundtrip_is_identity() {
        let v = 9;
        assert_eq!(Backend::from_mont(Backend::to_mont(v)), v);
    }

    #[test]
    fn beta_matches_configured_value() {
        assert_eq!(canonical(<Mod17Mont as Fp2Config>::BETA), BETA_CANONICAL);
    }

    #[test]
    fn const_to_mont_matches_runtime_to_mont() {
        // The compile-time path (to_mont_u32, used by Mod17Mont::BETA above)
        // must agree with the runtime path (MontWideBackend::to_mont) --
        // they implement the same REDC steps against the same inputs.
        assert_eq!(
            to_mont_u32(BETA_CANONICAL, Mod17Mont::R2, Mod17Mont::N_PRIME, Mod17Mont::MODULUS),
            Backend::to_mont(BETA_CANONICAL),
        );
    }

    #[test]
    fn add_matches_componentwise() {
        let sum = fe2(3, 5) + fe2(4, 9);
        assert_eq!(canonical(sum.c0), (3 + 4) % 17);
        assert_eq!(canonical(sum.c1), (5 + 9) % 17);
    }

    #[test]
    fn sub_is_inverse_of_add() {
        let a = fe2(3, 5);
        let b = fe2(4, 9);
        let back = (a + b) - b;
        assert_eq!(canonical(back.c0), canonical(a.c0));
        assert_eq!(canonical(back.c1), canonical(a.c1));
    }

    #[test]
    fn neg_matches_sub_from_zero() {
        let a = fe2(3, 5);
        let zero = fe2(0, 0);
        let neg = -a;
        assert_eq!(canonical(neg.c0), canonical((zero - a).c0));
        assert_eq!(canonical(neg.c1), canonical((zero - a).c1));
    }

    #[test]
    fn mul_matches_schoolbook_expansion() {
        let (a0, a1, b0, b1) = (3u32, 5u32, 4u32, 9u32);
        let product = fe2(a0, a1) * fe2(b0, b1);

        // (a0 + a1*u)(b0 + b1*u) = (a0*b0 + BETA*a1*b1) + (a0*b1 + a1*b0)*u
        let expected_c0 = (a0 * b0 + BETA_CANONICAL * a1 * b1) % 17;
        let expected_c1 = (a0 * b1 + a1 * b0) % 17;

        assert_eq!(canonical(product.c0), expected_c0);
        assert_eq!(canonical(product.c1), expected_c1);
    }

    #[test]
    fn norm_matches_definition() {
        let (c0, c1) = (3u32, 5u32);
        let norm = fe2(c0, c1).norm();

        // c0^2 - BETA*c1^2
        let p = 17u32;
        let expected = ((c0 * c0 + p * p - (BETA_CANONICAL * c1 * c1) % p) % p) % p;

        assert_eq!(canonical(norm), expected);
    }

    #[test]
    fn norm_matches_self_times_conjugate() {
        let a = fe2(3, 5);
        let conjugate = fe2(3, 0) - fe2(0, 5);
        let product = a * conjugate;

        // a * conjugate(a) lands purely in c0, equal to norm(a).
        assert_eq!(canonical(product.c0), canonical(a.norm()));
        assert_eq!(canonical(product.c1), 0);
    }

    proptest! {
        #[test]
        fn norm_matches_definition_proptest(c0 in 0u32..17, c1 in 0u32..17) {
            let norm = fe2(c0, c1).norm();

            let p = 17u32;
            let expected = ((c0 * c0 + p * p - (BETA_CANONICAL * c1 * c1) % p) % p) % p;

            prop_assert_eq!(canonical(norm), expected);
        }

        #[test]
        fn add_commutative(a0 in 0u32..17, a1 in 0u32..17, b0 in 0u32..17, b1 in 0u32..17) {
            let lhs = fe2(a0, a1) + fe2(b0, b1);
            let rhs = fe2(b0, b1) + fe2(a0, a1);
            prop_assert_eq!(canonical(lhs.c0), canonical(rhs.c0));
            prop_assert_eq!(canonical(lhs.c1), canonical(rhs.c1));
        }

        #[test]
        fn mul_commutative(a0 in 0u32..17, a1 in 0u32..17, b0 in 0u32..17, b1 in 0u32..17) {
            let lhs = fe2(a0, a1) * fe2(b0, b1);
            let rhs = fe2(b0, b1) * fe2(a0, a1);
            prop_assert_eq!(canonical(lhs.c0), canonical(rhs.c0));
            prop_assert_eq!(canonical(lhs.c1), canonical(rhs.c1));
        }

        #[test]
        fn mul_matches_schoolbook(a0 in 0u32..17, a1 in 0u32..17, b0 in 0u32..17, b1 in 0u32..17) {
            let product = fe2(a0, a1) * fe2(b0, b1);

            let expected_c0 = (a0 * b0 + BETA_CANONICAL * a1 * b1) % 17;
            let expected_c1 = (a0 * b1 + a1 * b0) % 17;

            prop_assert_eq!(canonical(product.c0), expected_c0);
            prop_assert_eq!(canonical(product.c1), expected_c1);
        }

        #[test]
        fn mul_distributes_over_add(a0 in 0u32..17, a1 in 0u32..17, b0 in 0u32..17, b1 in 0u32..17, c0 in 0u32..17, c1 in 0u32..17) {
            let a = fe2(a0, a1);
            let b = fe2(b0, b1);
            let c = fe2(c0, c1);

            let lhs = a * (b + c);
            let rhs = a * b + a * c;

            prop_assert_eq!(canonical(lhs.c0), canonical(rhs.c0));
            prop_assert_eq!(canonical(lhs.c1), canonical(rhs.c1));
        }

        #[test]
        fn inverse_times_self_is_one(c0 in 0u32..17, c1 in 0u32..17) {
            // BETA is a non-residue, so norm(a) == 0 only when a itself is
            // zero -- every other element is invertible.
            prop_assume!(c0 != 0 || c1 != 0);

            let a = fe2(c0, c1);
            let product = a * a.inverse();

            prop_assert_eq!(canonical(product.c0), 1);
            prop_assert_eq!(canonical(product.c1), 0);
        }
    }

    #[test]
    fn frobenius_coeff_matches_euler_criterion() {
        // BETA is a quadratic non-residue mod 17, so BETA^((p-1)/2) == -1.
        assert_eq!(canonical(<Mod17Mont as Fp2Config>::FROBENIUS_COEFF), 16);
    }

    #[test]
    fn frobenius_matches_conjugate_for_a_quadratic_extension() {
        // FROBENIUS_COEFF == -1 here (see the test above), so frobenius
        // collapses to conjugation -- true for any genuine Fp2, since
        // Euler's criterion always makes a non-residue's ((p-1)/2)th power
        // -1, never some other coefficient.
        let a = fe2(3, 5);
        let frob = a.frobenius();
        let conjugate = fe2(3, 0) - fe2(0, 5);

        assert_eq!(canonical(frob.c0), canonical(conjugate.c0));
        assert_eq!(canonical(frob.c1), canonical(conjugate.c1));
    }

    #[test]
    fn frobenius_fixes_the_base_field_embedding() {
        // c1 = 0 is Fp embedded in Fp2; Frobenius must fix it pointwise.
        let a = fe2(5, 0);
        let frob = a.frobenius();

        assert_eq!(canonical(frob.c0), canonical(a.c0));
        assert_eq!(canonical(frob.c1), 0);
    }

    proptest! {
        #[test]
        fn frobenius_matches_formula(c0 in 0u32..17, c1 in 0u32..17) {
            let frob = fe2(c0, c1).frobenius();

            let expected_c1 = (c1 * 16) % 17; // c1 * FROBENIUS_COEFF
            prop_assert_eq!(canonical(frob.c0), c0);
            prop_assert_eq!(canonical(frob.c1), expected_c1);
        }

        #[test]
        fn frobenius_is_an_involution(c0 in 0u32..17, c1 in 0u32..17) {
            // [Fp2:Fp] = 2, so applying Frobenius twice computes x^(p^2),
            // which is the identity on Fp2.
            let a = fe2(c0, c1);
            let twice = a.frobenius().frobenius();

            prop_assert_eq!(canonical(twice.c0), canonical(a.c0));
            prop_assert_eq!(canonical(twice.c1), canonical(a.c1));
        }

        #[test]
        fn frobenius_is_additive(a0 in 0u32..17, a1 in 0u32..17, b0 in 0u32..17, b1 in 0u32..17) {
            let a = fe2(a0, a1);
            let b = fe2(b0, b1);

            let lhs = (a + b).frobenius();
            let rhs = a.frobenius() + b.frobenius();

            prop_assert_eq!(canonical(lhs.c0), canonical(rhs.c0));
            prop_assert_eq!(canonical(lhs.c1), canonical(rhs.c1));
        }

        #[test]
        fn frobenius_is_multiplicative(a0 in 0u32..17, a1 in 0u32..17, b0 in 0u32..17, b1 in 0u32..17) {
            let a = fe2(a0, a1);
            let b = fe2(b0, b1);

            let lhs = (a * b).frobenius();
            let rhs = a.frobenius() * b.frobenius();

            prop_assert_eq!(canonical(lhs.c0), canonical(rhs.c0));
            prop_assert_eq!(canonical(lhs.c1), canonical(rhs.c1));
        }
    }

    // TODO: mont_add/sub/mul/neg-style cross-checks against a
    // non-Montgomery Fp2, mirroring field's DefaultBackend/MontBackend
    // parity tests.
}
