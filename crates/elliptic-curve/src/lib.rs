use std::fmt;
use std::ops::{Add, Mul};

use field_core::{Field, FpRepr};

// Short Weierstrass curve y^2 = x^3 + A*x + B over Field, with Scalar the
// (prime) order-n field used for scalar multiplication exponents -- kept
// distinct from Field since the two are almost always different fields
// (e.g. secp256k1's base field vs. its scalar field).
pub trait Curve {
    type Field: Field;
    type Scalar: Field;

    const A: Self::Field;
    const B: Self::Field;
}

// x, y generic over any Curve so the same point type works over Fp, Fp2,
// or the small fields alike. infinity marks the point at infinity (the
// group identity); x/y are meaningless when it's set, matching the usual
// affine short-Weierstrass convention.
pub struct AffinePoint<C: Curve> {
    pub x: C::Field,
    pub y: C::Field,
    pub infinity: bool,
}

// Derived impls would require C: Clone/Copy, but only C::Field needs to be
// (Field's Copy supertrait already guarantees that), so implement by hand --
// same reasoning as Fp<B>'s manual Clone/Copy in field-core.
impl<C: Curve> Clone for AffinePoint<C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<C: Curve> Copy for AffinePoint<C> {}

impl<C: Curve> fmt::Debug for AffinePoint<C>
where
    C::Field: fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AffinePoint")
            .field("x", &self.x)
            .field("y", &self.y)
            .field("infinity", &self.infinity)
            .finish()
    }
}

// No `where C::Field: PartialEq` needed: Field carries PartialEq as a
// supertrait, so C::Field: Field already guarantees it.
impl<C: Curve> PartialEq for AffinePoint<C> {
    fn eq(&self, other: &Self) -> bool {
        match (self.infinity, other.infinity) {
            (true, true) => true,
            (false, false) => self.x == other.x && self.y == other.y,
            _ => false,
        }
    }
}

// Standard short-Weierstrass chord-and-tangent addition. Leans on Field's
// PartialEq supertrait to tell the doubling case (self == other) apart
// from the vertical-line case (self == -other) when x1 == x2 -- both only
// differ by their y coordinate, so there's no way to route between the two
// formulas without comparing field elements.
impl<C: Curve> Add for AffinePoint<C> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        if self.infinity {
            return other;
        }
        if other.infinity {
            return self;
        }

        let lambda = if self.x == other.x {
            if self.y == -other.y {
                // self == -other: the chord is vertical, sum is infinity.
                // x/y are meaningless once infinity is set, so any values
                // do -- reuse self's, since we have them on hand.
                return Self { x: self.x, y: self.y, infinity: true };
            }

            // self == other: tangent slope, lambda = (3*x^2 + A) / (2*y).
            let x_sq = self.x.square();
            let numerator = x_sq + x_sq + x_sq + C::A;
            let denominator = self.y + self.y;

            numerator * denominator.inverse()
        } else {
            // Chord slope, lambda = (y2 - y1) / (x2 - x1).
            let numerator = other.y - self.y;
            let denominator = other.x - self.x;

            numerator * denominator.inverse()
        };

        let x3 = lambda.square() - self.x - other.x;
        let y3 = lambda * (self.x - x3) - self.y;

        Self { x: x3, y: y3, infinity: false }
    }
}

// Scalar multiplication by double-and-add: the additive-group analogue of
// Fp::pow's square-and-multiply (field-core/src/lib.rs) -- "double" takes
// the place of "square" and point "add" takes the place of field "mul",
// walking the scalar's bits from least to most significant. O(log scalar)
// point additions instead of O(scalar).
//
// The scalar is any R: FpRepr (a bit-indexable integer, e.g. u64 or
// bigint::U256) rather than C::Scalar itself: Field doesn't expose bit
// access on its elements (same reason Fp::pow's exponent is B::Repr, not
// Self), so there's no way to walk a C::Scalar's bits directly.
impl<C: Curve, R: FpRepr> Mul<R> for AffinePoint<C> {
    type Output = Self;

    fn mul(self, scalar: R) -> Self::Output {
        // Identity element O; x/y are meaningless once infinity is set,
        // so any values do -- reuse self's, same convention as Add.
        let mut result = Self { x: self.x, y: self.y, infinity: true };
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

    // Textbook toy curve y^2 = x^3 + 2x + 2 mod 17: order 19, generator
    // G = (5,1). Worked example values below (2G, 3G) are the standard
    // ones used to sanity-check short-Weierstrass addition by hand.
    struct Mod17;
    impl FpConfig for Mod17 {
        const MODULUS: U256 = U256::from_u64(17);
    }
    type F17 = Fp<DefaultBackend<Mod17>>;

    fn fe(v: u64) -> F17 {
        F17::new(U256::from(v % 17))
    }

    struct Curve17;
    impl Curve for Curve17 {
        type Field = F17;
        type Scalar = F17;

        const A: F17 = F17::new(U256::from_u64(2));
        const B: F17 = F17::new(U256::from_u64(2));
    }
    type P17 = AffinePoint<Curve17>;

    fn pt(x: u64, y: u64) -> P17 {
        P17 { x: fe(x), y: fe(y), infinity: false }
    }

    fn inf() -> P17 {
        P17 { x: fe(0), y: fe(0), infinity: true }
    }

    fn g() -> P17 {
        pt(5, 1)
    }

    #[test]
    fn identity_is_additive_identity() {
        assert!(g() + inf() == g());
        assert!(inf() + g() == g());
    }

    #[test]
    fn doubling_matches_worked_example() {
        assert!(g() + g() == pt(6, 3));
    }

    #[test]
    fn chord_addition_matches_worked_example() {
        let g2 = g() + g();
        assert!(g2 + g() == pt(10, 6));
    }

    #[test]
    fn point_plus_its_negation_is_infinity() {
        let neg_g = pt(5, 16); // -1 mod 17 == 16
        assert!(g() + neg_g == inf());
    }

    #[test]
    fn addition_is_commutative() {
        let g3 = g() + g() + g();
        assert!(g() + g3 == g3 + g());
    }

    #[test]
    fn generator_has_order_19() {
        let mut acc = inf();
        for _ in 0..19 {
            acc = acc + g();
        }
        assert!(acc == inf());
    }

    #[test]
    fn scalar_mul_by_zero_is_infinity() {
        assert!(g() * 0u64 == inf());
    }

    #[test]
    fn scalar_mul_by_one_is_self() {
        assert!(g() * 1u64 == g());
    }

    #[test]
    fn scalar_mul_matches_worked_examples() {
        assert!(g() * 2u64 == pt(6, 3));
        assert!(g() * 3u64 == pt(10, 6));
    }

    #[test]
    fn scalar_mul_matches_repeated_addition() {
        let mut acc = inf();
        for k in 0..25u64 {
            assert!(g() * k == acc, "mismatch at k={k}");
            acc = acc + g();
        }
    }

    #[test]
    fn scalar_mul_by_order_is_infinity() {
        assert!(g() * 19u64 == inf());
    }

    #[test]
    fn scalar_mul_wraps_modulo_order() {
        // 19 is the group order, so k and k+19 must land on the same point.
        assert!(g() * 5u64 == g() * 24u64);
    }

    #[test]
    fn scalar_mul_accepts_any_fprepr_not_just_u64() {
        // Same double-and-add path, driven by a bigint scalar instead of a
        // machine int -- confirms genericity over R: FpRepr (matters for
        // real curves, whose scalar field doesn't fit in a u64).
        assert!(g() * U256::from_u64(3) == pt(10, 6));
    }

    #[test]
    fn scalar_mul_of_infinity_is_infinity() {
        assert!(inf() * 7u64 == inf());
    }
}
