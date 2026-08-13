use std::fmt;
use std::ops::{Add, Mul};

use field_core::{Field, FpRepr};

// Short Weierstrass curve y^2 = x^3 + A*x + B over Field, with Scalar the
// (prime) order-n field used for scalar multiplication exponents -- kept
// distinct from Field since the two are almost always different fields
// (e.g. secp256k1's base field vs. its scalar field).
pub trait ShortWeierstrassCurve {
    type Field: Field;
    type Scalar: Field;

    const A: Self::Field;
    const B: Self::Field;
}

// x, y generic over any ShortWeierstrassCurve so the same point type works
// over Fp, Fp2, or the small fields alike. infinity marks the point at
// infinity (the group identity); x/y are meaningless when it's set,
// matching the usual affine short-Weierstrass convention.
pub struct AffinePoint<C: ShortWeierstrassCurve> {
    pub x: C::Field,
    pub y: C::Field,
    pub infinity: bool,
}

// Derived impls would require C: Clone/Copy, but only C::Field needs to be
// (Field's Copy supertrait already guarantees that), so implement by hand --
// same reasoning as Fp<B>'s manual Clone/Copy in field-core.
impl<C: ShortWeierstrassCurve> Clone for AffinePoint<C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<C: ShortWeierstrassCurve> Copy for AffinePoint<C> {}

impl<C: ShortWeierstrassCurve> fmt::Debug for AffinePoint<C>
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
impl<C: ShortWeierstrassCurve> PartialEq for AffinePoint<C> {
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
impl<C: ShortWeierstrassCurve> Add for AffinePoint<C> {
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

impl<C: ShortWeierstrassCurve> AffinePoint<C> {
    // Checks the defining equation y^2 = x^3 + A*x + B. The point at
    // infinity is the group's identity, not an affine solution to the
    // equation, so it's exempted and always validates.
    pub fn validate(&self) -> bool {
        if self.infinity {
            return true;
        }

        let lhs = self.y.square();
        let rhs = self.x.square() * self.x + C::A * self.x + C::B;

        lhs == rhs
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
//
// This is the cheaper of the two scalar_mul strategies (elliptic-curve's
// Cargo.toml `scalar-mul-*` features select exactly one at compile time),
// but both the number of additions and the branch taken each bit depend on
// the scalar -- see the constant-time ladder below for the alternative.
#[cfg(feature = "scalar-mul-variable-time")]
impl<C: ShortWeierstrassCurve, R: FpRepr> Mul<R> for AffinePoint<C> {
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

// Montgomery ladder, constant-time: maintains two accumulators r0 = k'*P,
// r1 = (k'+1)*P (k' being the bits of the scalar processed so far, from the
// top down), each iteration performing exactly one add and one double
// regardless of the bit -- unlike double-and-add above, which skips the
// add on a zero bit. Iterates a fixed `R::BITS` positions regardless of the
// scalar's magnitude (leading zero bits leave r0 = O, r1 = P unchanged, so
// this is still correct -- just walks extra no-op iterations for small
// scalars), and swaps r0/r1 via `cswap_affine` (field arithmetic, see
// ladder.rs) instead of a branch on the bit.
#[cfg(feature = "scalar-mul-constant-time")]
impl<C: ShortWeierstrassCurve, R: FpRepr> Mul<R> for AffinePoint<C> {
    type Output = Self;

    fn mul(self, scalar: R) -> Self::Output {
        let mut r0 = Self { x: self.x, y: self.y, infinity: true };
        let mut r1 = self;

        for i in (0..R::BITS).rev() {
            let bit = scalar.bit(i);

            let (a0, a1) = cswap_affine(bit, r0, r1);
            let sum = a0 + a1;
            let dbl = a0 + a0;
            let (b0, b1) = cswap_affine(bit, dbl, sum);
            r0 = b0;
            r1 = b1;
        }

        r0
    }
}

// Swaps (a, b) when bit is true, leaves them as-is when bit is false --
// via `ladder::select`/`select_bool` on each field instead of a branch on
// `bit`, which is the constant-time ladder's whole point (see its doc
// comment above).
#[cfg(feature = "scalar-mul-constant-time")]
fn cswap_affine<C: ShortWeierstrassCurve>(
    bit: bool,
    a: AffinePoint<C>,
    b: AffinePoint<C>,
) -> (AffinePoint<C>, AffinePoint<C>) {
    let out_a = AffinePoint {
        x: crate::ladder::select(bit, a.x, b.x),
        y: crate::ladder::select(bit, a.y, b.y),
        infinity: crate::ladder::select_bool(bit, a.infinity, b.infinity),
    };
    let out_b = AffinePoint {
        x: crate::ladder::select(bit, b.x, a.x),
        y: crate::ladder::select(bit, b.y, a.y),
        infinity: crate::ladder::select_bool(bit, b.infinity, a.infinity),
    };

    (out_a, out_b)
}

// Jacobian projective coordinates: (X, Y, Z) represents the affine point
// (X/Z^2, Y/Z^3). The point at infinity is any triple with Z == 0 (mirrors
// AffinePoint's `infinity` flag, but folded into the coordinate itself
// rather than a separate bool). The payoff over AffinePoint is that add/
// double need zero field inversions -- AffinePoint::add above needs one
// inverse per call -- at the cost of extra multiplications and a
// non-unique representation (the same affine point has many (X,Y,Z)
// triples), which is why PartialEq below cross-multiplies instead of
// comparing coordinates directly.
pub struct JacobianPoint<C: ShortWeierstrassCurve> {
    pub x: C::Field,
    pub y: C::Field,
    pub z: C::Field,
}

// Derived impls would require C: Clone/Copy, but only C::Field needs to be
// (Field's Copy supertrait already guarantees that), so implement by hand --
// same reasoning as AffinePoint's manual Clone/Copy above.
impl<C: ShortWeierstrassCurve> Clone for JacobianPoint<C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<C: ShortWeierstrassCurve> Copy for JacobianPoint<C> {}

impl<C: ShortWeierstrassCurve> fmt::Debug for JacobianPoint<C>
where
    C::Field: fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("JacobianPoint").field("x", &self.x).field("y", &self.y).field("z", &self.z).finish()
    }
}

// Jacobian coordinates aren't unique -- (X, Y, Z) and (r^2*X, r^3*Y, r*Z)
// denote the same affine point for any nonzero r -- so equality can't just
// compare coordinates like AffinePoint does. Instead cross-multiply: X1/Z1^2
// == X2/Z2^2 iff X1*Z2^2 == X2*Z1^2, same idea for Y, avoiding a field
// inversion just to test equality.
impl<C: ShortWeierstrassCurve> PartialEq for JacobianPoint<C> {
    fn eq(&self, other: &Self) -> bool {
        match (self.is_infinity(), other.is_infinity()) {
            (true, true) => true,
            (false, false) => {
                let z1z1 = self.z.square();
                let z2z2 = other.z.square();

                self.x * z2z2 == other.x * z1z1 && self.y * z2z2 * other.z == other.y * z1z1 * self.z
            }
            _ => false,
        }
    }
}

impl<C: ShortWeierstrassCurve> JacobianPoint<C> {
    // The point at infinity: any Z == 0 triple works (X/Z^2, Y/Z^3 are
    // meaningless once infinity is set, same convention as AffinePoint), so
    // (1, 1, 0) is as good as any other and matches the identity a fresh
    // "no point yet" accumulator would want.
    pub fn infinity() -> Self {
        Self { x: C::Field::one(), y: C::Field::one(), z: C::Field::zero() }
    }

    pub fn is_infinity(&self) -> bool {
        self.z == C::Field::zero()
    }

    pub fn from_affine(p: AffinePoint<C>) -> Self {
        if p.infinity {
            return Self::infinity();
        }

        Self { x: p.x, y: p.y, z: C::Field::one() }
    }

    // The one inversion Jacobian coordinates defer: recovers (x, y) =
    // (X/Z^2, Y/Z^3) from (X, Y, Z).
    pub fn to_affine(self) -> AffinePoint<C> {
        if self.is_infinity() {
            return AffinePoint { x: C::Field::zero(), y: C::Field::zero(), infinity: true };
        }

        let z_inv = self.z.inverse();
        let z_inv_sq = z_inv.square();

        AffinePoint { x: self.x * z_inv_sq, y: self.y * z_inv_sq * z_inv, infinity: false }
    }

    // Checks Y^2 = X^3 + A*X*Z^4 + B*Z^6 -- the short-Weierstrass equation
    // y^2 = x^3 + A*x + B with x = X/Z^2, y = Y/Z^3 cleared of denominators
    // by multiplying through by Z^6.
    pub fn validate(&self) -> bool {
        if self.is_infinity() {
            return true;
        }

        let zz = self.z.square();
        let zzzz = zz.square();

        let lhs = self.y.square();
        let rhs = self.x.square() * self.x + C::A * self.x * zzzz + C::B * zzzz * zz;

        lhs == rhs
    }

    // Point doubling ("dbl-2007-bl"): 4M + 6S, generic over the curve's A
    // coefficient (no a=0 / a=-3 shortcut) since ShortWeierstrassCurve
    // doesn't restrict A to those special-cased values.
    pub fn double(self) -> Self {
        if self.is_infinity() {
            return self;
        }

        let xx = self.x.square();
        let yy = self.y.square();
        let yyyy = yy.square();
        let zz = self.z.square();

        // s = 4*X*Y^2, via (X+YY)^2 - XX - YYYY = 2*X*YY (see the identity
        // this expands to: X^2 + 2*X*YY + YY^2 - XX - YYYY).
        let s = (self.x + yy).square() - xx - yyyy;
        let s = s + s;
        let m = xx + xx + xx + C::A * zz.square();
        let t = m.square() - s - s;

        let x3 = t;
        let yyyy8 = {
            let yyyy2 = yyyy + yyyy;
            let yyyy4 = yyyy2 + yyyy2;
            yyyy4 + yyyy4
        };
        let y3 = m * (s - t) - yyyy8;
        let z3 = (self.y + self.z).square() - yy - zz;

        Self { x: x3, y: y3, z: z3 }
    }
}

// Standard short-Weierstrass chord-and-tangent addition in Jacobian form
// ("add-2007-bl"), 12M + 4S plus one doubling on the coincident-point path.
// H == 0 marks the two points sharing an affine x-coordinate: r == 0 then
// means they're equal (routed to the cheaper dedicated `double`), otherwise
// they're mutual negatives (chord is vertical, sum is infinity) -- the same
// case split AffinePoint::add makes via direct x/y comparison, just
// detected through the cross-multiplied H/r terms instead since Jacobian
// coordinates aren't unique.
impl<C: ShortWeierstrassCurve> Add for JacobianPoint<C> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        if self.is_infinity() {
            return other;
        }
        if other.is_infinity() {
            return self;
        }

        let z1z1 = self.z.square();
        let z2z2 = other.z.square();
        let u1 = self.x * z2z2;
        let u2 = other.x * z1z1;
        let s1 = self.y * other.z * z2z2;
        let s2 = other.y * self.z * z1z1;

        let h = u2 - u1;
        let r = s2 - s1;

        if h == C::Field::zero() {
            if r == C::Field::zero() {
                return self.double();
            }
            return Self::infinity();
        }

        let i = (h + h).square();
        let j = h * i;
        let r = r + r;
        let v = u1 * i;

        let x3 = r.square() - j - (v + v);
        let y3 = r * (v - x3) - (s1 * j + s1 * j);
        let z3 = ((self.z + other.z).square() - z1z1 - z2z2) * h;

        Self { x: x3, y: y3, z: z3 }
    }
}

// Scalar multiplication by double-and-add, same structure as AffinePoint's
// Mul<R> impl above, but routed through the dedicated `double` (4M+6S)
// rather than `self + self` (12M+4S) for the squaring step -- the whole
// point of carrying Jacobian coordinates through scalar_mul instead of just
// add. Same `scalar-mul-*` feature selection as AffinePoint's Mul<R> above.
#[cfg(feature = "scalar-mul-variable-time")]
impl<C: ShortWeierstrassCurve, R: FpRepr> Mul<R> for JacobianPoint<C> {
    type Output = Self;

    fn mul(self, scalar: R) -> Self::Output {
        let mut result = Self::infinity();
        let mut base = self;
        let mut e = scalar;

        while e != R::ZERO {
            if e.bit(0) {
                result = result + base;
            }
            base = base.double();
            e >>= 1;
        }

        result
    }
}

// Montgomery ladder, constant-time -- same fixed-width `R::BITS` loop and
// arithmetic `cswap` (see ladder.rs) as AffinePoint's constant-time ladder
// above, routed through the dedicated `double` (4M+6S) instead of
// `self + self` (12M+4S), same payoff as the double-and-add impl above.
#[cfg(feature = "scalar-mul-constant-time")]
impl<C: ShortWeierstrassCurve, R: FpRepr> Mul<R> for JacobianPoint<C> {
    type Output = Self;

    fn mul(self, scalar: R) -> Self::Output {
        let mut r0 = Self::infinity();
        let mut r1 = self;

        for i in (0..R::BITS).rev() {
            let bit = scalar.bit(i);

            let (a0, a1) = cswap_jacobian(bit, r0, r1);
            let sum = a0 + a1;
            let dbl = a0.double();
            let (b0, b1) = cswap_jacobian(bit, dbl, sum);
            r0 = b0;
            r1 = b1;
        }

        r0
    }
}

// Swaps (a, b) when bit is true via `ladder::select` on each coordinate --
// no `infinity` bool to select here unlike AffinePoint's cswap, since
// JacobianPoint folds infinity into Z == 0 rather than a separate flag.
#[cfg(feature = "scalar-mul-constant-time")]
fn cswap_jacobian<C: ShortWeierstrassCurve>(
    bit: bool,
    a: JacobianPoint<C>,
    b: JacobianPoint<C>,
) -> (JacobianPoint<C>, JacobianPoint<C>) {
    let out_a = JacobianPoint {
        x: crate::ladder::select(bit, a.x, b.x),
        y: crate::ladder::select(bit, a.y, b.y),
        z: crate::ladder::select(bit, a.z, b.z),
    };
    let out_b = JacobianPoint {
        x: crate::ladder::select(bit, b.x, a.x),
        y: crate::ladder::select(bit, b.y, a.y),
        z: crate::ladder::select(bit, b.z, a.z),
    };

    (out_a, out_b)
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
    impl ShortWeierstrassCurve for Curve17 {
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

    #[test]
    fn validate_accepts_points_on_the_curve() {
        // G, 2G, 3G from the worked example, plus its negation -- all
        // genuine solutions to y^2 = x^3 + 2x + 2 mod 17.
        assert!(g().validate());
        assert!(pt(6, 3).validate());
        assert!(pt(10, 6).validate());
        assert!(pt(5, 16).validate());
    }

    #[test]
    fn validate_rejects_points_off_the_curve() {
        // Same x as G but the wrong y: 5^3 + 2*5 + 2 = 1 mod 17, whose
        // square roots are 1 and 16, not 2.
        assert!(!pt(5, 2).validate());
        // x^3+2x+2 mod 17 at x=0 is 2 (whose square roots are 6 and 11,
        // not 0), so (0,0) isn't on the curve either.
        assert!(!pt(0, 0).validate());
    }

    #[test]
    fn validate_accepts_infinity() {
        assert!(inf().validate());
    }

    #[test]
    fn validate_accepts_every_point_produced_by_add_and_scalar_mul() {
        let mut acc = inf();
        for k in 0..25u64 {
            assert!(acc.validate(), "k={k}");
            assert!((g() * k).validate(), "k={k}");
            acc = acc + g();
        }
    }

    type J17 = JacobianPoint<Curve17>;

    fn jg() -> J17 {
        J17::from_affine(g())
    }

    fn jinf() -> J17 {
        J17::infinity()
    }

    // Rescales (X, Y, Z) by an arbitrary nonzero factor r -- a different
    // Jacobian representative of the same affine point, used to check that
    // PartialEq/to_affine are representation-independent rather than
    // comparing raw coordinates.
    fn rescale(p: J17, r: u64) -> J17 {
        let r = fe(r);
        let r2 = r.square();
        let r3 = r2 * r;
        J17 { x: p.x * r2, y: p.y * r3, z: p.z * r }
    }

    #[test]
    fn jacobian_infinity_round_trips_through_affine() {
        assert!(jinf().to_affine() == inf());
        assert!(J17::from_affine(inf()).is_infinity());
    }

    #[test]
    fn jacobian_affine_round_trip_preserves_point() {
        assert!(jg().to_affine() == g());
    }

    #[test]
    fn jacobian_rescaled_representations_are_equal() {
        assert!(jg() == rescale(jg(), 3));
        assert!(rescale(jg(), 3) == rescale(jg(), 11));
    }

    #[test]
    fn jacobian_identity_is_additive_identity() {
        assert!(jg() + jinf() == jg());
        assert!(jinf() + jg() == jg());
    }

    #[test]
    fn jacobian_doubling_matches_affine_doubling() {
        let expected = g() + g();
        assert!((jg() + jg()).to_affine() == expected);
        assert!(jg().double().to_affine() == expected);
    }

    #[test]
    fn jacobian_doubling_agrees_with_generic_add() {
        // The dedicated `double` (4M+6S) must land on the same point as
        // routing self+self through the general chord-and-tangent formula.
        assert!(jg().double() == jg() + jg());
    }

    #[test]
    fn jacobian_chord_addition_matches_affine_addition() {
        let g2j = jg() + jg();
        let g2a = g() + g();
        assert!((g2j + jg()).to_affine() == g2a + g());
    }

    #[test]
    fn jacobian_point_plus_its_negation_is_infinity() {
        // AffinePoint's own negation, -(x, y) = (x, -y) -- Field::Neg is a
        // supertrait, so this needs no hardcoded coordinate.
        let neg_g = J17::from_affine(AffinePoint { x: g().x, y: -g().y, infinity: false });
        assert!(jg() + neg_g == jinf());
    }

    #[test]
    fn jacobian_add_matches_affine_add_across_many_pairs() {
        // Sweeps k, m over the generator's full order (19), so this covers
        // doubling (k == m), mutual negation (k + m == 19), and either
        // operand being infinity (k == 0 or m == 0) against AffinePoint's
        // independent Add implementation -- not just the two worked
        // examples above.
        for k in 0..19u64 {
            for m in 0..19u64 {
                let expected = (g() * k) + (g() * m);
                let actual = ((jg() * k) + (jg() * m)).to_affine();
                assert!(actual == expected, "mismatch at k={k}, m={m}");
            }
        }
    }

    #[test]
    fn jacobian_addition_is_commutative() {
        let g3 = jg() + jg() + jg();
        assert!(jg() + g3 == g3 + jg());
    }

    #[test]
    fn jacobian_addition_matches_rescaled_operands() {
        let g2 = jg() + jg();
        assert!(rescale(jg(), 7) + rescale(g2, 13) == jg() + g2);
    }

    #[test]
    fn jacobian_generator_has_order_19() {
        let mut acc = jinf();
        for _ in 0..19 {
            acc = acc + jg();
        }
        assert!(acc == jinf());
    }

    #[test]
    fn jacobian_scalar_mul_by_zero_is_infinity() {
        assert!(jg() * 0u64 == jinf());
    }

    #[test]
    fn jacobian_scalar_mul_matches_affine_scalar_mul() {
        // Cross-check against AffinePoint's independent scalar_mul
        // implementation over the same curve and scalar, including scalars
        // past the generator's order (19) to exercise the wraparound too.
        for k in 0..25u64 {
            assert!((jg() * k).to_affine() == g() * k, "mismatch at k={k}");
        }
    }

    #[test]
    fn jacobian_scalar_mul_accepts_any_fprepr_not_just_u64() {
        assert!(jg() * U256::from_u64(3) == jg() * 3u64);
    }

    #[test]
    fn jacobian_validate_accepts_points_on_the_curve() {
        assert!(jg().validate());
        assert!((jg() * 2u64).validate());
        assert!((jg() * 3u64).validate());
        assert!(rescale(jg(), 7).validate());
    }

    #[test]
    fn jacobian_validate_rejects_points_off_the_curve() {
        assert!(!J17::from_affine(pt(5, 2)).validate());
    }

    #[test]
    fn jacobian_validate_accepts_infinity() {
        assert!(jinf().validate());
    }

    #[test]
    fn jacobian_validate_accepts_every_point_produced_by_add_and_scalar_mul() {
        let mut acc = jinf();
        for k in 0..25u64 {
            assert!(acc.validate(), "k={k}");
            assert!((jg() * k).validate(), "k={k}");
            acc = acc + jg();
        }
    }
}
