use std::fmt;
use std::ops::Add;

use field_core::Field;

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

impl<C: Curve> PartialEq for AffinePoint<C>
where
    C::Field: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        match (self.infinity, other.infinity) {
            (true, true) => true,
            (false, false) => self.x == other.x && self.y == other.y,
            _ => false,
        }
    }
}

// Standard short-Weierstrass chord-and-tangent addition. Requires
// C::Field: PartialEq to tell the doubling case (self == other) apart from
// the vertical-line case (self == -other) when x1 == x2 -- both only
// differ by their y coordinate, so there's no way to route between the two
// formulas without comparing field elements.
impl<C: Curve> Add for AffinePoint<C>
where
    C::Field: PartialEq,
{
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
