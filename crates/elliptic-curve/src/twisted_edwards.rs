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
// place of the point-at-infinity accumulator start. Selected at compile
// time by elliptic-curve's Cargo.toml `scalar-mul-*` features, same as
// every other point type's Mul<R> impl in this crate -- see
// short_weierstrass.rs's AffinePoint Mul<R> impls for the full rationale
// shared across both variants below.
#[cfg(feature = "scalar-mul-variable-time")]
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

// Montgomery ladder, constant-time -- see
// short_weierstrass::AffinePoint's constant-time ladder impl for the full
// rationale; identical shape here, just with Edwards' identity() in place
// of the infinity-flagged accumulator start.
#[cfg(feature = "scalar-mul-constant-time")]
impl<C: TwistedEdwardsCurve, R: FpRepr> Mul<R> for AffinePoint<C> {
    type Output = Self;

    fn mul(self, scalar: R) -> Self::Output {
        let mut r0 = Self::identity();
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

// Swaps (a, b) when bit is true via `ladder::select` on each coordinate --
// no bool flag to select here, since twisted Edwards' identity is a genuine
// affine point rather than an out-of-band infinity marker (see
// AffinePoint's doc comment above).
#[cfg(feature = "scalar-mul-constant-time")]
fn cswap_affine<C: TwistedEdwardsCurve>(
    bit: bool,
    a: AffinePoint<C>,
    b: AffinePoint<C>,
) -> (AffinePoint<C>, AffinePoint<C>) {
    let out_a = AffinePoint { x: crate::ladder::select(bit, a.x, b.x), y: crate::ladder::select(bit, a.y, b.y) };
    let out_b = AffinePoint { x: crate::ladder::select(bit, b.x, a.x), y: crate::ladder::select(bit, b.y, a.y) };

    (out_a, out_b)
}

// Extended twisted Edwards coordinates (Hisil-Wong-Carter-Dawson 2008):
// (X, Y, Z, T) represents the affine point (X/Z, Y/Z), with the redundant T
// coordinate held equal to X*Y/Z (equivalently X*Y == Z*T) so the unified
// addition law below never needs to recompute the x1*y1*x2*y2-style cross
// terms AffinePoint::add's D::D * x1*x2*y1*y2 needs. Mirrors
// short_weierstrass::JacobianPoint's payoff over AffinePoint: add/double
// need zero field inversions -- AffinePoint::add above needs two -- at the
// cost of extra multiplications and a non-unique representation (the same
// affine point has many (X,Y,Z,T) quadruples), which is why PartialEq below
// cross-multiplies instead of comparing coordinates directly. Unlike
// JacobianPoint there's no separate infinity flag or Z==0 sentinel to
// special-case: the identity (0,1) lifts to (0, Z, Z, 0) for any nonzero Z,
// which the same formulas that handle every other point already produce
// and consume correctly (same "no exceptional cases" completeness
// AffinePoint::add relies on).
pub struct ExtendedPoint<C: TwistedEdwardsCurve> {
    pub x: C::Field,
    pub y: C::Field,
    pub z: C::Field,
    pub t: C::Field,
}

// Derived impls would require C: Clone/Copy, but only C::Field needs to be
// (Field's Copy supertrait already guarantees that), so implement by hand --
// same reasoning as AffinePoint's manual Clone/Copy above.
impl<C: TwistedEdwardsCurve> Clone for ExtendedPoint<C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<C: TwistedEdwardsCurve> Copy for ExtendedPoint<C> {}

impl<C: TwistedEdwardsCurve> fmt::Debug for ExtendedPoint<C>
where
    C::Field: fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ExtendedPoint")
            .field("x", &self.x)
            .field("y", &self.y)
            .field("z", &self.z)
            .field("t", &self.t)
            .finish()
    }
}

// Extended coordinates aren't unique -- (X, Y, Z, T) and (r*X, r*Y, r*Z,
// r*T) denote the same affine point for any nonzero r -- so equality can't
// just compare coordinates like AffinePoint does. Instead cross-multiply:
// X1/Z1 == X2/Z2 iff X1*Z2 == X2*Z1, same idea for Y, avoiding a field
// inversion just to test equality (same approach as JacobianPoint::eq,
// simpler here since extended coordinates scale x/y linearly rather than by
// Z^2/Z^3).
impl<C: TwistedEdwardsCurve> PartialEq for ExtendedPoint<C> {
    fn eq(&self, other: &Self) -> bool {
        self.x * other.z == other.x * self.z && self.y * other.z == other.y * self.z
    }
}

impl<C: TwistedEdwardsCurve> ExtendedPoint<C> {
    // The group identity, (0, 1) lifted to extended coordinates: (0, 1, 1, 0)
    // satisfies X*Y == Z*T (0*1 == 1*0) and matches AffinePoint::identity
    // under the X/Z, Y/Z projection.
    pub fn identity() -> Self {
        Self { x: C::Field::zero(), y: C::Field::one(), z: C::Field::one(), t: C::Field::zero() }
    }

    pub fn from_affine(p: AffinePoint<C>) -> Self {
        Self { x: p.x, y: p.y, z: C::Field::one(), t: p.x * p.y }
    }

    // The one inversion extended coordinates defer: recovers (x, y) =
    // (X/Z, Y/Z) from (X, Y, Z, T).
    pub fn to_affine(self) -> AffinePoint<C> {
        let z_inv = self.z.inverse();

        AffinePoint { x: self.x * z_inv, y: self.y * z_inv }
    }

    // Checks A*X^2 + Y^2 = Z^2 + D*T^2 -- the twisted Edwards equation
    // A*x^2 + y^2 = 1 + D*x^2*y^2 with x = X/Z, y = Y/Z cleared of
    // denominators by multiplying through by Z^2, then x^2*y^2 = (X*Y/Z)^2
    // rewritten as T^2 via the extended-coordinate invariant X*Y = Z*T. That
    // invariant is itself checked separately, since a quadruple could
    // otherwise satisfy the curve equation with a T that doesn't actually
    // track X*Y/Z.
    pub fn validate(&self) -> bool {
        if self.x * self.y != self.z * self.t {
            return false;
        }

        let x_sq = self.x.square();
        let y_sq = self.y.square();
        let z_sq = self.z.square();
        let t_sq = self.t.square();

        C::A * x_sq + y_sq == z_sq + C::D * t_sq
    }

    // Dedicated doubling ("dbl-2008-hwcd"), generic over the curve's A
    // coefficient: 4M + 4S + 1 mul-by-A, cheaper than routing self + self
    // through the general unified `Add` below (which spends extra
    // multiplications computing cross terms that collapse when both inputs
    // are the same point).
    pub fn double(self) -> Self {
        let xx = self.x.square();
        let yy = self.y.square();
        let zz2 = self.z.square() + self.z.square();
        let d = C::A * xx;
        let e = (self.x + self.y).square() - xx - yy;
        let g = d + yy;
        let f = g - zz2;
        let h = d - yy;

        Self { x: e * f, y: g * h, z: f * g, t: e * h }
    }
}

// Unified extended-coordinate addition ("add-2008-hwcd-4"), generic over the
// curve's A coefficient: 9M + 1 mul-by-D + 1 mul-by-A. "Unified" the same
// way AffinePoint::add above is -- this single formula also handles
// doubling (self == other) and the identity on either side, with no
// exceptional cases to branch on, for the same complete-curve reason
// AffinePoint::add's doc comment gives (A a nonzero square, D a
// non-square). Unlike JacobianPoint::add there's no H==0/r==0 case split
// needed before calling `double`: completeness means this formula alone is
// always correct, so `double` above exists purely as a cheaper shortcut,
// not a correctness requirement.
impl<C: TwistedEdwardsCurve> Add for ExtendedPoint<C> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        let a = self.x * other.x;
        let b = self.y * other.y;
        let c = C::D * self.t * other.t;
        let d = self.z * other.z;
        let e = (self.x + self.y) * (other.x + other.y) - a - b;
        let f = d - c;
        let g = d + c;
        let h = b - C::A * a;

        Self { x: e * f, y: g * h, z: f * g, t: e * h }
    }
}

// Scalar multiplication by double-and-add, same structure as AffinePoint's
// Mul<R> impl, but routed through the dedicated `double` (4M+4S) rather
// than `self + self` (9M) for the squaring step -- the whole point of
// carrying extended coordinates through scalar_mul instead of just add,
// mirroring JacobianPoint::mul's reasoning. Same `scalar-mul-*` feature
// selection as every other Mul<R> impl in this crate.
#[cfg(feature = "scalar-mul-variable-time")]
impl<C: TwistedEdwardsCurve, R: FpRepr> Mul<R> for ExtendedPoint<C> {
    type Output = Self;

    fn mul(self, scalar: R) -> Self::Output {
        let mut result = Self::identity();
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

// Montgomery ladder, constant-time -- see
// short_weierstrass::JacobianPoint's constant-time ladder impl for the full
// rationale; identical shape here, routed through the dedicated `double`
// for the same reason the double-and-add impl above is.
#[cfg(feature = "scalar-mul-constant-time")]
impl<C: TwistedEdwardsCurve, R: FpRepr> Mul<R> for ExtendedPoint<C> {
    type Output = Self;

    fn mul(self, scalar: R) -> Self::Output {
        let mut r0 = Self::identity();
        let mut r1 = self;

        for i in (0..R::BITS).rev() {
            let bit = scalar.bit(i);

            let (a0, a1) = cswap_extended(bit, r0, r1);
            let sum = a0 + a1;
            let dbl = a0.double();
            let (b0, b1) = cswap_extended(bit, dbl, sum);
            r0 = b0;
            r1 = b1;
        }

        r0
    }
}

// Swaps (a, b) when bit is true via `ladder::select` on each of the four
// coordinates -- no bool flag to select, same reason as
// twisted_edwards::AffinePoint's cswap above.
#[cfg(feature = "scalar-mul-constant-time")]
fn cswap_extended<C: TwistedEdwardsCurve>(
    bit: bool,
    a: ExtendedPoint<C>,
    b: ExtendedPoint<C>,
) -> (ExtendedPoint<C>, ExtendedPoint<C>) {
    let out_a = ExtendedPoint {
        x: crate::ladder::select(bit, a.x, b.x),
        y: crate::ladder::select(bit, a.y, b.y),
        z: crate::ladder::select(bit, a.z, b.z),
        t: crate::ladder::select(bit, a.t, b.t),
    };
    let out_b = ExtendedPoint {
        x: crate::ladder::select(bit, b.x, a.x),
        y: crate::ladder::select(bit, b.y, a.y),
        z: crate::ladder::select(bit, b.z, a.z),
        t: crate::ladder::select(bit, b.t, a.t),
    };

    (out_a, out_b)
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

    type E41 = ExtendedPoint<Curve41>;

    fn eg() -> E41 {
        E41::from_affine(g())
    }

    fn eidentity() -> E41 {
        E41::identity()
    }

    // Rescales (X, Y, Z, T) by an arbitrary nonzero factor r -- a different
    // extended representative of the same affine point, used to check that
    // PartialEq/to_affine are representation-independent rather than
    // comparing raw coordinates.
    fn rescale(p: E41, r: i64) -> E41 {
        let r = fe(r);
        E41 { x: p.x * r, y: p.y * r, z: p.z * r, t: p.t * r }
    }

    #[test]
    fn extended_identity_round_trips_through_affine() {
        assert!(eidentity().to_affine() == P41::identity());
        assert!(E41::from_affine(P41::identity()) == eidentity());
    }

    #[test]
    fn extended_affine_round_trip_preserves_point() {
        assert!(eg().to_affine() == g());
    }

    #[test]
    fn extended_rescaled_representations_are_equal() {
        assert!(eg() == rescale(eg(), 3));
        assert!(rescale(eg(), 3) == rescale(eg(), 11));
    }

    #[test]
    fn extended_identity_is_additive_identity() {
        assert!(eg() + eidentity() == eg());
        assert!(eidentity() + eg() == eg());
    }

    #[test]
    fn extended_doubling_matches_affine_doubling() {
        let expected = g() + g();
        assert!((eg() + eg()).to_affine() == expected);
        assert!(eg().double().to_affine() == expected);
    }

    #[test]
    fn extended_doubling_agrees_with_generic_add() {
        // The dedicated `double` (4M+4S) must land on the same point as
        // routing self+self through the general unified addition formula.
        assert!(eg().double() == eg() + eg());
    }

    #[test]
    fn extended_chord_addition_matches_affine_addition() {
        let g2e = eg() + eg();
        let g2a = g() + g();
        assert!((g2e + eg()).to_affine() == g2a + g());
    }

    #[test]
    fn extended_point_plus_its_negation_is_identity() {
        // Edwards negation: -(x, y) = (-x, y).
        let neg_g = E41::from_affine(pt(-4, 2));
        assert!(eg() + neg_g == eidentity());
    }

    #[test]
    fn extended_addition_is_commutative() {
        let g3 = eg() + eg() + eg();
        assert!(eg() + g3 == g3 + eg());
    }

    #[test]
    fn extended_addition_matches_rescaled_operands() {
        let g2 = eg() + eg();
        assert!(rescale(eg(), 7) + rescale(g2, 13) == eg() + g2);
    }

    #[test]
    fn extended_scalar_mul_by_zero_is_identity() {
        assert!(eg() * 0u64 == eidentity());
    }

    #[test]
    fn extended_scalar_mul_by_one_is_self() {
        assert!(eg() * 1u64 == eg());
    }

    #[test]
    fn extended_scalar_mul_matches_repeated_addition() {
        let mut acc = eidentity();
        for k in 0..15u64 {
            assert!(eg() * k == acc, "mismatch at k={k}");
            acc = acc + eg();
        }
    }

    #[test]
    fn extended_scalar_mul_accepts_any_fprepr_not_just_u64() {
        assert!(eg() * U256::from_u64(3) == eg() * 3u64);
    }

    #[test]
    fn extended_scalar_mul_matches_affine_scalar_mul() {
        // Cross-check against AffinePoint's independent scalar_mul
        // implementation over the same curve, across the generator's full
        // order (48, per the Curve41 fixture comment above).
        for k in 0..48u64 {
            assert!((eg() * k).to_affine() == g() * k, "mismatch at k={k}");
        }
    }

    #[test]
    fn extended_add_matches_affine_add_across_many_pairs() {
        // Sweeps k, m over the generator's full order (48), covering
        // doubling (k == m), mutual negation (k + m == 48), and either
        // operand being the identity (k == 0 or m == 0) against
        // AffinePoint's independent Add implementation.
        for k in 0..48u64 {
            for m in 0..48u64 {
                let expected = (g() * k) + (g() * m);
                let actual = ((eg() * k) + (eg() * m)).to_affine();
                assert!(actual == expected, "mismatch at k={k}, m={m}");
            }
        }
    }

    #[test]
    fn extended_generator_has_order_48() {
        let mut acc = eidentity();
        for _ in 0..48 {
            acc = acc + eg();
        }
        assert!(acc == eidentity());
    }

    #[test]
    fn extended_validate_accepts_points_on_the_curve() {
        assert!(eg().validate());
        assert!((eg() * 2u64).validate());
        assert!((eg() * 3u64).validate());
        assert!(rescale(eg(), 7).validate());
    }

    #[test]
    fn extended_validate_accepts_identity() {
        assert!(eidentity().validate());
    }

    #[test]
    fn extended_validate_rejects_points_off_the_curve() {
        assert!(!E41::from_affine(pt(1, 1)).validate());
    }

    #[test]
    fn extended_validate_rejects_inconsistent_t() {
        // Same X, Y, Z as a genuine curve point, but T doesn't track X*Y/Z.
        let mut p = eg();
        p.t = p.t + F41::one();
        assert!(!p.validate());
    }

    #[test]
    fn extended_validate_accepts_every_point_produced_by_add_and_scalar_mul() {
        let mut acc = eidentity();
        for k in 0..15u64 {
            assert!(acc.validate(), "k={k}");
            assert!((eg() * k).validate(), "k={k}");
            acc = acc + eg();
        }
    }
}
