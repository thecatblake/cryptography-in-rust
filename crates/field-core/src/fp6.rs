use std::ops::{Add, Mul, Neg, Sub};

use crate::{Fp2, Fp2Config};

// Fp6Config extends Fp2Config with the one extra constant a cubic extension
// needs: Fp6 = Fp2[v] / (v^3 - XI). Fp2Config already supplies everything
// needed to build the base field's Fp2 arithmetic, so an Fp6Config
// implementor only adds XI on top.
pub trait Fp6Config: Fp2Config + Sized {
    // Must be a cubic non-residue in Fp2 (XI has no cube root in Fp2), or
    // v^3 - XI factors and Fp6 collapses into Fp2 x Fp2 x Fp2 instead of
    // being a field.
    const XI: Fp2<Self>;
}

// Fp6 = Fp2[v] / (v^3 - XI), elements represented as c0 + c1*v + c2*v^2.
pub struct Fp6<C: Fp6Config> {
    pub c0: Fp2<C>,
    pub c1: Fp2<C>,
    pub c2: Fp2<C>,
}

// Derived impls would require C: Clone/Copy, but only C::Repr needs to be
// (FpRepr's Copy supertrait already guarantees that), so implement by hand
// -- same reasoning as Fp2<C>'s hand-written Clone/Copy.
impl<C: Fp6Config> Clone for Fp6<C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<C: Fp6Config> Copy for Fp6<C> {}

impl<C: Fp6Config> Fp6<C> {
    pub fn new(c0: Fp2<C>, c1: Fp2<C>, c2: Fp2<C>) -> Self {
        Fp6 { c0, c1, c2 }
    }

    // Multiplication by a fixed a = a0 + a1*v + a2*v^2 is the linear map on
    // (Fp2)^3 given, in the {1, v, v^2} basis and using v^3 = XI, by the
    // matrix (column j = a * basis_vector[j], reduced mod v^3 - XI)
    //   M = [ a0      XI*a2   XI*a1 ]
    //       [ a1      a0      XI*a2 ]
    //       [ a2      a1      a0    ]
    // a^-1 is the vector x solving M x = e1 = (1,0,0), i.e. x = M^-1 e1 =
    // adj(M) e1 / det(M) -- the first column of adj(M) (equivalently, M's
    // first-row cofactors) scaled by 1/det(M).
    //
    // Expanding det(M) along its first row and reusing the three minors
    // between the numerator and the determinant (worked out ahead of time
    // rather than run as a generic Cramer's-rule expansion of all nine 2x2
    // minors):
    //   t0 = a0^2 - XI*a1*a2   (cofactor of M[0][0], = adj(M) row 0 col 0)
    //   t1 = XI*a2^2 - a0*a1   (cofactor of M[0][1])
    //   t2 = a1^2 - a0*a2      (cofactor of M[0][2])
    //   det(M) = a0*t0 + XI*(a1*t2 + a2*t1)
    // and a^-1 = (t0, t1, t2) / det(M): one Fp2 inversion (of the scalar
    // det(M)) plus three Fp2 muls to scale, instead of inverting each of
    // the nine minors' worth of arithmetic separately. Total: 9 Fp2 muls
    // (6 to build t0..t2, 3 to fold them into det) + 3 XI const-muls + 3
    // Fp2 muls to scale by det^-1, plus the one Fp2 inversion.
    pub fn inverse(self) -> Self {
        let (a0, a1, a2) = (self.c0, self.c1, self.c2);

        let t0 = a0 * a0 - C::XI * (a1 * a2);
        let t1 = C::XI * (a2 * a2) - a0 * a1;
        let t2 = a1 * a1 - a0 * a2;

        let det = a0 * t0 + C::XI * (a1 * t2 + a2 * t1);
        let det_inv = det.inverse();

        Fp6::new(t0 * det_inv, t1 * det_inv, t2 * det_inv)
    }
}

impl<C: Fp6Config> Add for Fp6<C> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Fp6::new(self.c0 + rhs.c0, self.c1 + rhs.c1, self.c2 + rhs.c2)
    }
}

impl<C: Fp6Config> Sub for Fp6<C> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Fp6::new(self.c0 - rhs.c0, self.c1 - rhs.c1, self.c2 - rhs.c2)
    }
}

impl<C: Fp6Config> Neg for Fp6<C> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Fp6::new(-self.c0, -self.c1, -self.c2)
    }
}

// Karatsuba for a cubic extension: schoolbook
// (a0+a1*v+a2*v^2)(b0+b1*v+b2*v^2) needs 9 Fp2 muls. Its unreduced
// polynomial product has coefficients
//   d0 = a0*b0
//   d1 = a0*b1 + a1*b0
//   d2 = a0*b2 + a1*b1 + a2*b0
//   d3 = a1*b2 + a2*b1
//   d4 = a2*b2
// Naming v0=a0*b0, v1=a1*b1, v2=a2*b2, each cross-term pair collapses to one
// extra mul the same way Fp2::mul folds its own cross terms:
//   d1 = (a0+a1)*(b0+b1) - v0 - v1
//   d3 = (a1+a2)*(b1+b2) - v1 - v2
//   d2 = (a0+a2)*(b0+b2) - v0 - v2 + v1
// so d0..d4 cost 6 Fp2 muls total (v0, v1, v2, and the three cross
// products) instead of 9. Reducing mod v^3 - XI folds v^3 -> XI and
// v^4 -> XI*v:
//   c0 = d0 + XI*d3
//   c1 = d1 + XI*d4 = d1 + XI*v2
//   c2 = d2
// The two XI multiplies are against the fixed constant XI, not one of the
// two operands being multiplied together, so (as with Fp2::mul's BETA
// multiply) they're kept separate from the 6-mul count above.
impl<C: Fp6Config> Mul for Fp6<C> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let v0 = self.c0 * rhs.c0;
        let v1 = self.c1 * rhs.c1;
        let v2 = self.c2 * rhs.c2;

        let d1 = (self.c0 + self.c1) * (rhs.c0 + rhs.c1) - v0 - v1;
        let d3 = (self.c1 + self.c2) * (rhs.c1 + rhs.c2) - v1 - v2;
        let d2 = (self.c0 + self.c2) * (rhs.c0 + rhs.c2) - v0 - v2 + v1;

        let c0 = v0 + C::XI * d3;
        let c1 = d1 + C::XI * v2;
        let c2 = d2;

        Fp6::new(c0, c1, c2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Fp, MontFieldConfig, MontWideBackend, to_mont_u32};
    use proptest::prelude::*;

    // Same Mod17Mont fixture as fp2's tests: p = 17, R = 2^32.
    struct Mod17Mont;
    impl MontFieldConfig for Mod17Mont {
        type Repr = u32;

        const MODULUS: u32 = 17;
        const R2: u32 = 0x1;
        const N_PRIME: u32 = 0x0f0f_0f0f;
    }

    const BETA_CANONICAL: u32 = 3;

    impl Fp2Config for Mod17Mont {
        const BETA: Fp<MontWideBackend<Self>> =
            Fp::new(to_mont_u32(BETA_CANONICAL, Self::R2, Self::N_PRIME, Self::MODULUS));
    }

    // Canonical (non-Montgomery) coordinates of Mod17Mont::XI, used by tests
    // that check Fp6 arithmetic against a plain-integer formula.
    const XI_CANONICAL: (u32, u32) = (1, 1);

    impl Fp6Config for Mod17Mont {
        // Verified offline (see the accompanying derivation script): among
        // the 288 = 17^2-1 nonzero elements of F17^2, only 96 are cubes, and
        // 1+u isn't one of them, so v^3 - (1+u) is irreducible over F17^2.
        const XI: Fp2<Self> = Fp2 {
            c0: Fp::new(to_mont_u32(XI_CANONICAL.0, Self::R2, Self::N_PRIME, Self::MODULUS)),
            c1: Fp::new(to_mont_u32(XI_CANONICAL.1, Self::R2, Self::N_PRIME, Self::MODULUS)),
        };
    }

    type Backend = MontWideBackend<Mod17Mont>;
    type F17 = Fp<Backend>;
    type F17_2 = Fp2<Mod17Mont>;
    type F17_6 = Fp6<Mod17Mont>;

    fn fe(v: u32) -> F17 {
        F17::new(Backend::to_mont(v % 17))
    }

    fn canonical(x: F17) -> u32 {
        Backend::from_mont(x.value)
    }

    fn fe2(c0: u32, c1: u32) -> F17_2 {
        Fp2::new(fe(c0), fe(c1))
    }

    fn fe6(c0: (u32, u32), c1: (u32, u32), c2: (u32, u32)) -> F17_6 {
        Fp6::new(fe2(c0.0, c0.1), fe2(c1.0, c1.1), fe2(c2.0, c2.1))
    }

    fn canonical2(x: F17_2) -> (u32, u32) {
        (canonical(x.c0), canonical(x.c1))
    }

    #[test]
    fn xi_matches_configured_value() {
        assert_eq!(canonical2(Mod17Mont::XI), XI_CANONICAL);
    }

    #[test]
    fn add_matches_componentwise() {
        let sum = fe6((1, 2), (3, 4), (5, 6)) + fe6((7, 8), (9, 10), (11, 12));
        assert_eq!(canonical2(sum.c0), (8 % 17, 10 % 17));
        assert_eq!(canonical2(sum.c1), (12 % 17, 14 % 17));
        assert_eq!(canonical2(sum.c2), (16 % 17, 18 % 17));
    }

    #[test]
    fn sub_is_inverse_of_add() {
        let a = fe6((1, 2), (3, 4), (5, 6));
        let b = fe6((7, 8), (9, 10), (11, 12));
        let back = (a + b) - b;
        assert_eq!(canonical2(back.c0), canonical2(a.c0));
        assert_eq!(canonical2(back.c1), canonical2(a.c1));
        assert_eq!(canonical2(back.c2), canonical2(a.c2));
    }

    #[test]
    fn neg_matches_sub_from_zero() {
        let a = fe6((1, 2), (3, 4), (5, 6));
        let zero = fe6((0, 0), (0, 0), (0, 0));
        let neg = -a;
        assert_eq!(canonical2(neg.c0), canonical2((zero - a).c0));
        assert_eq!(canonical2(neg.c1), canonical2((zero - a).c1));
        assert_eq!(canonical2(neg.c2), canonical2((zero - a).c2));
    }

    // Plain-integer schoolbook reference for Fp2 multiplication, used only
    // by the Fp6 schoolbook check below.
    fn mul2_ref(a: (u32, u32), b: (u32, u32)) -> (u32, u32) {
        let p = 17u32;
        let c0 = (a.0 * b.0 + BETA_CANONICAL * a.1 * b.1) % p;
        let c1 = (a.0 * b.1 + a.1 * b.0) % p;
        (c0, c1)
    }

    fn add2_ref(a: (u32, u32), b: (u32, u32)) -> (u32, u32) {
        ((a.0 + b.0) % 17, (a.1 + b.1) % 17)
    }

    #[test]
    fn mul_matches_schoolbook_expansion() {
        let a = ((1u32, 2u32), (3u32, 4u32), (5u32, 6u32));
        let b = ((7u32, 8u32), (9u32, 10u32), (11u32, 12u32));

        let product = fe6(a.0, a.1, a.2) * fe6(b.0, b.1, b.2);

        // Unreduced polynomial product, then fold v^3 -> XI, v^4 -> XI*v.
        let d0 = mul2_ref(a.0, b.0);
        let d1 = add2_ref(mul2_ref(a.0, b.1), mul2_ref(a.1, b.0));
        let d2 = add2_ref(add2_ref(mul2_ref(a.0, b.2), mul2_ref(a.1, b.1)), mul2_ref(a.2, b.0));
        let d3 = add2_ref(mul2_ref(a.1, b.2), mul2_ref(a.2, b.1));
        let d4 = mul2_ref(a.2, b.2);

        let expected_c0 = add2_ref(d0, mul2_ref(XI_CANONICAL, d3));
        let expected_c1 = add2_ref(d1, mul2_ref(XI_CANONICAL, d4));
        let expected_c2 = d2;

        assert_eq!(canonical2(product.c0), expected_c0);
        assert_eq!(canonical2(product.c1), expected_c1);
        assert_eq!(canonical2(product.c2), expected_c2);
    }

    proptest! {
        #[test]
        fn add_commutative(
            a0 in 0u32..17, a1 in 0u32..17, a2 in 0u32..17, a3 in 0u32..17, a4 in 0u32..17, a5 in 0u32..17,
            b0 in 0u32..17, b1 in 0u32..17, b2 in 0u32..17, b3 in 0u32..17, b4 in 0u32..17, b5 in 0u32..17,
        ) {
            let a = fe6((a0, a1), (a2, a3), (a4, a5));
            let b = fe6((b0, b1), (b2, b3), (b4, b5));
            let lhs = a + b;
            let rhs = b + a;
            prop_assert_eq!(canonical2(lhs.c0), canonical2(rhs.c0));
            prop_assert_eq!(canonical2(lhs.c1), canonical2(rhs.c1));
            prop_assert_eq!(canonical2(lhs.c2), canonical2(rhs.c2));
        }

        #[test]
        fn mul_commutative(
            a0 in 0u32..17, a1 in 0u32..17, a2 in 0u32..17, a3 in 0u32..17, a4 in 0u32..17, a5 in 0u32..17,
            b0 in 0u32..17, b1 in 0u32..17, b2 in 0u32..17, b3 in 0u32..17, b4 in 0u32..17, b5 in 0u32..17,
        ) {
            let a = fe6((a0, a1), (a2, a3), (a4, a5));
            let b = fe6((b0, b1), (b2, b3), (b4, b5));
            let lhs = a * b;
            let rhs = b * a;
            prop_assert_eq!(canonical2(lhs.c0), canonical2(rhs.c0));
            prop_assert_eq!(canonical2(lhs.c1), canonical2(rhs.c1));
            prop_assert_eq!(canonical2(lhs.c2), canonical2(rhs.c2));
        }

        #[test]
        fn mul_distributes_over_add(
            a0 in 0u32..17, a1 in 0u32..17, a2 in 0u32..17, a3 in 0u32..17, a4 in 0u32..17, a5 in 0u32..17,
            b0 in 0u32..17, b1 in 0u32..17, b2 in 0u32..17, b3 in 0u32..17, b4 in 0u32..17, b5 in 0u32..17,
            c0 in 0u32..17, c1 in 0u32..17, c2 in 0u32..17, c3 in 0u32..17, c4 in 0u32..17, c5 in 0u32..17,
        ) {
            let a = fe6((a0, a1), (a2, a3), (a4, a5));
            let b = fe6((b0, b1), (b2, b3), (b4, b5));
            let c = fe6((c0, c1), (c2, c3), (c4, c5));

            let lhs = a * (b + c);
            let rhs = a * b + a * c;

            prop_assert_eq!(canonical2(lhs.c0), canonical2(rhs.c0));
            prop_assert_eq!(canonical2(lhs.c1), canonical2(rhs.c1));
            prop_assert_eq!(canonical2(lhs.c2), canonical2(rhs.c2));
        }

        #[test]
        fn inverse_times_self_is_one(
            a0 in 0u32..17, a1 in 0u32..17, a2 in 0u32..17, a3 in 0u32..17, a4 in 0u32..17, a5 in 0u32..17,
        ) {
            // Zero is the only element that can fail to be invertible; skip
            // it rather than asserting anything about Fp6::inverse there.
            prop_assume!(a0 != 0 || a1 != 0 || a2 != 0 || a3 != 0 || a4 != 0 || a5 != 0);

            let a = fe6((a0, a1), (a2, a3), (a4, a5));
            let product = a * a.inverse();

            prop_assert_eq!(canonical2(product.c0), (1, 0));
            prop_assert_eq!(canonical2(product.c1), (0, 0));
            prop_assert_eq!(canonical2(product.c2), (0, 0));
        }
    }
}
