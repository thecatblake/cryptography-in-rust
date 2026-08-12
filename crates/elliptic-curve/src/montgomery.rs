use std::fmt;
use std::ops::{Add, Mul};

use field_core::{Field, FpRepr};

// Montgomery curve B*y^2 = x^3 + A*x^2 + x over Field, with Scalar the
// order-n scalar field, same split as ShortWeierstrassCurve /
// TwistedEdwardsCurve. Kept as its own trait rather than reusing either:
// the defining equation and its addition/doubling formulas are shaped
// differently from both (a B coefficient scaling y^2 that neither of the
// other two forms has), so a point type for it needs its own arithmetic.
pub trait MontgomeryCurve {
    type Field: Field;
    type Scalar: Field;

    const A: Self::Field;
    const B: Self::Field;
}

// x, y generic over any MontgomeryCurve so the same point type works over
// Fp, Fp2, or the small fields alike. infinity marks the point at infinity
// (the group identity) -- same convention as short_weierstrass::AffinePoint,
// since a Montgomery curve's identity is likewise not a genuine affine
// solution to the defining equation; x/y are meaningless when it's set.
pub struct AffinePoint<C: MontgomeryCurve> {
    pub x: C::Field,
    pub y: C::Field,
    pub infinity: bool,
}

// Derived impls would require C: Clone/Copy, but only C::Field needs to be
// (Field's Copy supertrait already guarantees that), so implement by hand --
// same reasoning as short_weierstrass::AffinePoint's manual Clone/Copy.
impl<C: MontgomeryCurve> Clone for AffinePoint<C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<C: MontgomeryCurve> Copy for AffinePoint<C> {}

impl<C: MontgomeryCurve> fmt::Debug for AffinePoint<C>
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
impl<C: MontgomeryCurve> PartialEq for AffinePoint<C> {
    fn eq(&self, other: &Self) -> bool {
        match (self.infinity, other.infinity) {
            (true, true) => true,
            (false, false) => self.x == other.x && self.y == other.y,
            _ => false,
        }
    }
}

// Standard Montgomery chord-and-tangent addition. Same x1 == x2 branch
// structure as short_weierstrass::AffinePoint's Add: compares y to tell
// doubling (self == other) apart from the vertical-line case
// (self == -other), since both only differ by their y coordinate.
impl<C: MontgomeryCurve> Add for AffinePoint<C> {
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

            // self == other: tangent slope,
            // lambda = (3*x^2 + 2*A*x + 1) / (2*B*y).
            let x_sq = self.x.square();
            let numerator = x_sq + x_sq + x_sq + (C::A * self.x + C::A * self.x) + C::Field::one();
            let denominator = C::B * (self.y + self.y);

            numerator * denominator.inverse()
        } else {
            // Chord slope, lambda = (y2 - y1) / (x2 - x1).
            let numerator = other.y - self.y;
            let denominator = other.x - self.x;

            numerator * denominator.inverse()
        };

        let lambda_sq = lambda.square();
        let x3 = C::B * lambda_sq - C::A - self.x - other.x;
        let y3 = lambda * (self.x - x3) - self.y;

        Self { x: x3, y: y3, infinity: false }
    }
}

impl<C: MontgomeryCurve> AffinePoint<C> {
    // Checks the defining equation B*y^2 = x^3 + A*x^2 + x. The point at
    // infinity is the group's identity, not an affine solution to the
    // equation, so it's exempted and always validates.
    pub fn validate(&self) -> bool {
        if self.infinity {
            return true;
        }

        let lhs = C::B * self.y.square();
        let rhs = self.x.square() * self.x + C::A * self.x.square() + self.x;

        lhs == rhs
    }
}

// Scalar multiplication by double-and-add, same structure as
// short_weierstrass::AffinePoint's Mul<R> impl: "double" takes the place
// of Fp::pow's "square" and point "add" takes the place of field "mul",
// walking the scalar's bits from least to most significant.
impl<C: MontgomeryCurve, R: FpRepr> Mul<R> for AffinePoint<C> {
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

    // Toy curve y^2 = x^3 + 3x^2 + x mod 17 (A = 3, B = 1): order 16,
    // generator G = (5,1). Worked example values below (2G, 3G, 4G) were
    // computed independently by repeated application of the addition law.
    struct Mod17;
    impl FpConfig for Mod17 {
        const MODULUS: U256 = U256::from_u64(17);
    }
    type F17 = Fp<DefaultBackend<Mod17>>;

    fn fe(v: u64) -> F17 {
        F17::new(U256::from(v % 17))
    }

    struct Curve17;
    impl MontgomeryCurve for Curve17 {
        type Field = F17;
        type Scalar = F17;

        const A: F17 = F17::new(U256::from_u64(3));
        const B: F17 = F17::new(U256::from_u64(1));
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
        assert!(g() + g() == pt(8, 10));
    }

    #[test]
    fn chord_addition_matches_worked_example() {
        let g2 = g() + g();
        assert!(g2 + g() == pt(10, 1));
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
    fn generator_has_order_16() {
        let mut acc = inf();
        for _ in 0..16 {
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
        assert!(g() * 2u64 == pt(8, 10));
        assert!(g() * 3u64 == pt(10, 1));
        assert!(g() * 4u64 == pt(16, 16));
    }

    #[test]
    fn scalar_mul_matches_repeated_addition() {
        let mut acc = inf();
        for k in 0..20u64 {
            assert!(g() * k == acc, "mismatch at k={k}");
            acc = acc + g();
        }
    }

    #[test]
    fn scalar_mul_by_order_is_infinity() {
        assert!(g() * 16u64 == inf());
    }

    #[test]
    fn scalar_mul_wraps_modulo_order() {
        // 16 is the group order, so k and k+16 must land on the same point.
        assert!(g() * 5u64 == g() * 21u64);
    }

    #[test]
    fn scalar_mul_accepts_any_fprepr_not_just_u64() {
        // Same double-and-add path, driven by a bigint scalar instead of a
        // machine int -- confirms genericity over R: FpRepr (matters for
        // real curves, whose scalar field doesn't fit in a u64).
        assert!(g() * U256::from_u64(3) == pt(10, 1));
    }

    #[test]
    fn scalar_mul_of_infinity_is_infinity() {
        assert!(inf() * 7u64 == inf());
    }

    #[test]
    fn validate_accepts_points_on_the_curve() {
        // G, 2G, 3G, 4G from the worked example, plus its negation -- all
        // genuine solutions to y^2 = x^3 + 3x^2 + x mod 17.
        assert!(g().validate());
        assert!(pt(8, 10).validate());
        assert!(pt(10, 1).validate());
        assert!(pt(16, 16).validate());
        assert!(pt(5, 16).validate());
    }

    #[test]
    fn validate_rejects_points_off_the_curve() {
        // Same x as G but the wrong y: 5^3 + 3*5^2 + 5 = 1 mod 17, whose
        // square roots are 1 and 16, not 2.
        assert!(!pt(5, 2).validate());
    }

    #[test]
    fn validate_accepts_infinity() {
        assert!(inf().validate());
    }

    #[test]
    fn validate_accepts_every_point_produced_by_add_and_scalar_mul() {
        let mut acc = inf();
        for k in 0..20u64 {
            assert!(acc.validate(), "k={k}");
            assert!((g() * k).validate(), "k={k}");
            acc = acc + g();
        }
    }
}
