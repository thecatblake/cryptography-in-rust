use field_core::Field;

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
