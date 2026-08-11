use std::fmt;

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
