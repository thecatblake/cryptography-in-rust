use std::fmt;
use std::ops::{Add, Mul};

use field_core::{Field, FpRepr};

// Twisted Edwards curve A*x^2 + y^2 = 1 + D*x^2*y^2 over Field, with Scalar
// the order-n scalar field, same split as ShortWeierstrassCurve. Distinct
// trait rather than a shared "two curve constants" trait: the constants
// mean different things (A/B are the Weierstrass cubic's coefficients; A/D
// shape a totally different quartic-in-disguise curve with its own
// addition law and identity point at (0,1) rather than at infinity), so a
// point type for it would need its own arithmetic rather than reusing
// short_weierstrass::AffinePoint's chord-and-tangent formulas.
pub trait TwistedEdwardsCurve {
    type Field: Field;
    type Scalar: Field;

    const A: Self::Field;
    const D: Self::Field;
}

// x, y generic over any TwistedEdwardsCurve so the same point type works
// over Fp, Fp2, or the small fields alike. Unlike short_weierstrass's
// AffinePoint, there's no separate infinity flag: the identity is the
// genuine affine point (0, 1), which the unified addition law below
// produces and consumes like any other point.
pub struct AffinePoint<C: TwistedEdwardsCurve> {
    pub x: C::Field,
    pub y: C::Field,
}

// Derived impls would require C: Clone/Copy, but only C::Field needs to be
// (Field's Copy supertrait already guarantees that), so implement by hand --
// same reasoning as short_weierstrass::AffinePoint's manual Clone/Copy.
impl<C: TwistedEdwardsCurve> Clone for AffinePoint<C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<C: TwistedEdwardsCurve> Copy for AffinePoint<C> {}

impl<C: TwistedEdwardsCurve> fmt::Debug for AffinePoint<C>
where
    C::Field: fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AffinePoint").field("x", &self.x).field("y", &self.y).finish()
    }
}

// No `where C::Field: PartialEq` needed: Field carries PartialEq as a
// supertrait, so C::Field: Field already guarantees it.
impl<C: TwistedEdwardsCurve> PartialEq for AffinePoint<C> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<C: TwistedEdwardsCurve> AffinePoint<C> {
    // The group identity, (0, 1): substituting into A*x^2 + y^2 =
    // 1 + D*x^2*y^2 gives 0 + 1 = 1 + 0, true for any A, D.
    pub fn identity() -> Self {
        Self { x: C::Field::zero(), y: C::Field::one() }
    }

    // Checks the defining equation A*x^2 + y^2 = 1 + D*x^2*y^2.
    pub fn validate(&self) -> bool {
        let x_sq = self.x.square();
        let y_sq = self.y.square();

        let lhs = C::A * x_sq + y_sq;
        let rhs = C::Field::one() + C::D * x_sq * y_sq;

        lhs == rhs
    }
}

// Unified twisted Edwards addition:
//   x3 = (x1*y2 + y1*x2) / (1 + D*x1*x2*y1*y2)
//   y3 = (y1*y2 - A*x1*x2) / (1 - D*x1*x2*y1*y2)
// "Unified" means this single formula also handles doubling (self == other)
// -- unlike short Weierstrass, there's no separate tangent-line case to
// branch on. This holds for every pair of points precisely when the curve
// is complete (D is a non-square in Field and A is a nonzero square),
// which guarantees both denominators are always nonzero; TwistedEdwardsCurve
// doesn't enforce that algebraically; it's on the curve's chosen A/D to
// satisfy it.
impl<C: TwistedEdwardsCurve> Add for AffinePoint<C> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        let one = C::Field::one();

        let x1y2 = self.x * other.y;
        let y1x2 = self.y * other.x;
        let y1y2 = self.y * other.y;
        let x1x2 = self.x * other.x;
        let d_prod = C::D * x1x2 * y1y2;

        let x3 = (x1y2 + y1x2) * (one + d_prod).inverse();
        let y3 = (y1y2 - C::A * x1x2) * (one - d_prod).inverse();

        Self { x: x3, y: y3 }
    }
}

// Scalar multiplication by double-and-add, same structure as
// short_weierstrass::AffinePoint's Mul<R> impl: "double" (self-addition)
// takes the place of Fp::pow's "square", and Edwards' identity() takes the
// place of the point-at-infinity accumulator start.
impl<C: TwistedEdwardsCurve, R: FpRepr> Mul<R> for AffinePoint<C> {
    type Output = Self;

    fn mul(self, scalar: R) -> Self::Output {
        let mut result = Self::identity();
        let mut base = self;
        let mut e = scalar;

        while e != R::ZERO {
            if e.bit(0) {
                result = result + base;
            }
            base = base + base;
            e >>= 1;
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bigint::U256;
    use field::{DefaultBackend, Fp, FpConfig};

    // Toy curve -x^2 + y^2 = 1 - 3*x^2*y^2 mod 41 (A = -1 = 40, D = 3).
    // 41 == 1 mod 4, so -1 is a quadratic residue mod 41 (satisfying A being
    // a nonzero square), and 3 was checked (by brute-force search over
    // 1..41) to be a quadratic non-residue mod 41 (satisfying D being a
    // non-square) -- so the curve is complete: every pair of affine points
    // adds via the single unified formula with no exceptional cases.
    // (4, 2) was found the same way, checked against the curve equation and
    // confirmed (by repeated addition) to generate the full 48-point group.
    struct Mod41;
    impl FpConfig for Mod41 {
        const MODULUS: U256 = U256::from_u64(41);
    }
    type F41 = Fp<DefaultBackend<Mod41>>;

    fn fe(v: i64) -> F41 {
        F41::new(U256::from(v.rem_euclid(41) as u64))
    }

    struct Curve41;
    impl TwistedEdwardsCurve for Curve41 {
        type Field = F41;
        type Scalar = F41;

        const A: F41 = F41::new(U256::from_u64(40));
        const D: F41 = F41::new(U256::from_u64(3));
    }
    type P41 = AffinePoint<Curve41>;

    fn pt(x: i64, y: i64) -> P41 {
        P41 { x: fe(x), y: fe(y) }
    }

    fn g() -> P41 {
        pt(4, 2)
    }

    #[test]
    fn generator_is_on_curve() {
        assert!(g().validate());
    }

    #[test]
    fn identity_is_on_curve() {
        assert!(P41::identity().validate());
    }

    #[test]
    fn identity_is_additive_identity() {
        assert!(g() + P41::identity() == g());
        assert!(P41::identity() + g() == g());
    }

    #[test]
    fn addition_is_commutative() {
        let g2 = g() + g();
        let g3 = g2 + g();
        assert!(g() + g3 == g3 + g());
    }

    #[test]
    fn doubling_matches_worked_example() {
        // Computed independently by repeated application of the addition
        // law starting from (0,1).
        assert!(g() + g() == pt(26, 19));
    }

    #[test]
    fn chord_addition_matches_worked_example() {
        let g2 = g() + g();
        assert!(g2 + g() == pt(16, 31));
    }

    #[test]
    fn point_plus_its_negation_is_identity() {
        // Edwards negation: -(x, y) = (-x, y).
        let neg_g = pt(-4, 2);
        assert!(g() + neg_g == P41::identity());
    }

    #[test]
    fn scalar_mul_by_zero_is_identity() {
        assert!(g() * 0u64 == P41::identity());
    }

    #[test]
    fn scalar_mul_by_one_is_self() {
        assert!(g() * 1u64 == g());
    }

    #[test]
    fn scalar_mul_matches_repeated_addition() {
        let mut acc = P41::identity();
        for k in 0..15u64 {
            assert!(g() * k == acc, "mismatch at k={k}");
            acc = acc + g();
        }
    }

    #[test]
    fn scalar_mul_accepts_any_fprepr_not_just_u64() {
        let expected = g() + g() + g();
        assert!(g() * U256::from_u64(3) == expected);
    }

    #[test]
    fn validate_accepts_every_point_produced_by_add_and_scalar_mul() {
        let mut acc = P41::identity();
        for k in 0..15u64 {
            assert!(acc.validate(), "k={k}");
            assert!((g() * k).validate(), "k={k}");
            acc = acc + g();
        }
    }

    #[test]
    fn validate_rejects_points_off_the_curve() {
        assert!(!pt(1, 1).validate());
    }
}
