use std::ops::{Add, Mul, Neg, Sub};

use crate::{Fp, MontFieldConfig, MontWideBackend};

// Fp2Config extends MontFieldConfig with the one extra constant a quadratic
// extension needs: Fp2 = Fp[u] / (u^2 - BETA). The MontFieldConfig
// supertrait already supplies everything needed to build the base field's
// backend (MontWideBackend<Self>, via MODULUS/R2/N_PRIME), so an Fp2Config
// implementor holds both the base field's Montgomery backend and BETA.
//
// BETA is stored in canonical (non-Montgomery) form, same as MODULUS --
// `beta()` below converts it once it's needed as a field element.
pub trait Fp2Config: MontFieldConfig {
    // Must be a quadratic non-residue mod Self::MODULUS, or u^2 - BETA
    // factors and Fp2 collapses to Fp x Fp instead of being a field.
    const BETA: Self::Repr;
}

// Fp2 = Fp[u] / (u^2 - BETA), elements represented as c0 + c1*u.
pub struct Fp2<C: Fp2Config> {
    pub c0: Fp<MontWideBackend<C>>,
    pub c1: Fp<MontWideBackend<C>>,
}

// Derived impls would require C: Clone/Copy, but only C::Repr needs to be
// (FpRepr's Copy supertrait already guarantees that), so implement by hand
// -- same reasoning as Fp<B>'s hand-written Clone/Copy in lib.rs.
impl<C: Fp2Config> Clone for Fp2<C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<C: Fp2Config> Copy for Fp2<C> {}

impl<C: Fp2Config> Fp2<C> {
    pub fn new(c0: Fp<MontWideBackend<C>>, c1: Fp<MontWideBackend<C>>) -> Self {
        Fp2 { c0, c1 }
    }

    // BETA lifted into the base field, in Montgomery form.
    pub fn beta() -> Fp<MontWideBackend<C>> {
        Fp::new(MontWideBackend::<C>::to_mont(C::BETA))
    }

    // N(a) = a * conjugate(a) = (c0+c1*u)(c0-c1*u) = c0^2 - c1^2*u^2, and
    // u^2 == BETA by definition of Fp2, so this collapses to a base-field
    // element. `inverse` below is defined in terms of it.
    pub fn norm(self) -> Fp<MontWideBackend<C>> {
        self.c0.square() - Self::beta() * self.c1.square()
    }

    // Multiplicative inverse. For a = c0 + c1*u,
    // a^-1 = conjugate(a) / norm(a) = (c0 - c1*u) * norm(a)^-1. The
    // denominator is the same c0^2 - BETA*c1^2 as `norm`, computed in place
    // here rather than via self.norm() since it's needed before the
    // conjugate is built.
    pub fn inverse(self) -> Self {
        let denom_inv = (self.c0.square() - Self::beta() * self.c1.square()).inverse();

        Fp2::new(self.c0 * denom_inv, -(self.c1 * denom_inv))
    }
}

impl<C: Fp2Config> Add for Fp2<C> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Fp2::new(self.c0 + rhs.c0, self.c1 + rhs.c1)
    }
}

impl<C: Fp2Config> Sub for Fp2<C> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Fp2::new(self.c0 - rhs.c0, self.c1 - rhs.c1)
    }
}

impl<C: Fp2Config> Neg for Fp2<C> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Fp2::new(-self.c0, -self.c1)
    }
}

// Karatsuba: schoolbook (a0+a1*u)(b0+b1*u) needs 4 base-field muls
// (a0*b0, a1*b1, a0*b1, a1*b0). Naming v0 = a0*b0 and v1 = a1*b1, the two
// cross terms' sum a0*b1 + a1*b0 equals (a0+a1)*(b0+b1) - v0 - v1, so only
// one more mul -- not two -- is needed to get both of them at once:
//   c0 = v0 + BETA*v1
//   c1 = (a0+a1)*(b0+b1) - v0 - v1
// 3 base-field muls total (v0, v1, and the cross-term product); the BETA
// multiply is separate since BETA is a fixed constant, not one of the two
// operands being multiplied together.
impl<C: Fp2Config> Mul for Fp2<C> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let v0 = self.c0 * rhs.c0;
        let v1 = self.c1 * rhs.c1;

        let c0 = v0 + Self::beta() * v1;
        let c1 = (self.c0 + self.c1) * (rhs.c0 + rhs.c1) - v0 - v1;

        Fp2::new(c0, c1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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

    impl Fp2Config for Mod17Mont {
        // QRs mod 17 are {1,2,4,8,9,13,15,16}; 3 isn't among them, so
        // u^2 - 3 is irreducible over F17.
        const BETA: u32 = 3;
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
        assert_eq!(canonical(F17_2::beta()), Mod17Mont::BETA);
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
        let expected_c0 = (a0 * b0 + Mod17Mont::BETA * a1 * b1) % 17;
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
        let expected = ((c0 * c0 + p * p - (Mod17Mont::BETA * c1 * c1) % p) % p) % p;

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
            let expected = ((c0 * c0 + p * p - (Mod17Mont::BETA * c1 * c1) % p) % p) % p;

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

            let expected_c0 = (a0 * b0 + Mod17Mont::BETA * a1 * b1) % 17;
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

    // TODO: mont_add/sub/mul/neg-style cross-checks against a
    // non-Montgomery Fp2, mirroring field's DefaultBackend/MontBackend
    // parity tests.
}
