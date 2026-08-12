use std::ops::{Add, Mul, Neg, Sub};

use crate::{Field, Frobenius};

// CubicExtConfig describes a cubic extension Base[v] / (v^3 - XI) for some
// base field Base. Fp6Config (Base = Fp2<C>, i.e. QuadExt<C>) is the only
// instantiation used elsewhere in this crate, but as with QuadExtConfig,
// Base only needs to be a Field, not Fp2 specifically.
pub trait CubicExtConfig: Sized {
    type Base: Field;

    // Must be a cubic non-residue in Base (XI has no cube root in Base), or
    // v^3 - XI factors and the extension collapses into Base x Base x Base
    // instead of being a field.
    const XI: Self::Base;
}

// CubicExtFrobeniusConfig supplies the two extra pieces of data Frobenius
// needs beyond CubicExtConfig: XI^((p-1)/3) and XI^(2(p-1)/3) in Base, so
// that v^p and v^(2p) reduce to FROBENIUS_COEFF_C1*v and
// FROBENIUS_COEFF_C2*v^2 respectively (see the Frobenius impl below).
// Split out from CubicExtConfig itself for the same reason
// QuadExtFrobeniusConfig is split from QuadExtConfig: computing these
// coefficients needs the base field's actual characteristic, and in
// practice compile-time exponentiation machinery, which not every
// CubicExtConfig implementor should be forced to provide. Fp6Config is
// this crate's only implementor (fp6.rs blanket-impls this trait for any
// C: Fp6Config).
pub trait CubicExtFrobeniusConfig: CubicExtConfig {
    const FROBENIUS_COEFF_C1: Self::Base;
    const FROBENIUS_COEFF_C2: Self::Base;
}

// CubicExt = Base[v] / (v^3 - XI), elements represented as c0 + c1*v + c2*v^2.
pub struct CubicExt<C: CubicExtConfig> {
    pub c0: C::Base,
    pub c1: C::Base,
    pub c2: C::Base,
}

// Derived impls would require C: Clone/Copy, but only C::Base needs to be
// (Field's Copy supertrait already guarantees that), so implement by hand
// -- same reasoning as QuadExt<C>'s hand-written Clone/Copy.
impl<C: CubicExtConfig> Clone for CubicExt<C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<C: CubicExtConfig> Copy for CubicExt<C> {}

// No `where C::Base: PartialEq` needed: C::Base: Field, and Field now
// carries PartialEq as a supertrait.
impl<C: CubicExtConfig> PartialEq for CubicExt<C> {
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1 && self.c2 == other.c2
    }
}

impl<C: CubicExtConfig> CubicExt<C> {
    pub fn new(c0: C::Base, c1: C::Base, c2: C::Base) -> Self {
        CubicExt { c0, c1, c2 }
    }

    // Multiplication by a fixed a = a0 + a1*v + a2*v^2 is the linear map on
    // Base^3 given, in the {1, v, v^2} basis and using v^3 = XI, by the
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
    // and a^-1 = (t0, t1, t2) / det(M): one Base inversion (of the scalar
    // det(M)) plus three Base muls to scale, instead of inverting each of
    // the nine minors' worth of arithmetic separately. Total: 9 Base muls
    // (6 to build t0..t2, 3 to fold them into det) + 3 XI const-muls + 3
    // Base muls to scale by det^-1, plus the one Base inversion.
    pub fn inverse(self) -> Self {
        let (a0, a1, a2) = (self.c0, self.c1, self.c2);

        let t0 = a0.square() - C::XI * (a1 * a2);
        let t1 = C::XI * a2.square() - a0 * a1;
        let t2 = a1.square() - a0 * a2;

        let det = a0 * t0 + C::XI * (a1 * t2 + a2 * t1);
        let det_inv = det.inverse();

        CubicExt::new(t0 * det_inv, t1 * det_inv, t2 * det_inv)
    }
}

impl<C: CubicExtConfig> Add for CubicExt<C> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        CubicExt::new(self.c0 + rhs.c0, self.c1 + rhs.c1, self.c2 + rhs.c2)
    }
}

impl<C: CubicExtConfig> Sub for CubicExt<C> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        CubicExt::new(self.c0 - rhs.c0, self.c1 - rhs.c1, self.c2 - rhs.c2)
    }
}

impl<C: CubicExtConfig> Neg for CubicExt<C> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        CubicExt::new(-self.c0, -self.c1, -self.c2)
    }
}

// Karatsuba for a cubic extension: schoolbook
// (a0+a1*v+a2*v^2)(b0+b1*v+b2*v^2) needs 9 Base muls. Its unreduced
// polynomial product has coefficients
//   d0 = a0*b0
//   d1 = a0*b1 + a1*b0
//   d2 = a0*b2 + a1*b1 + a2*b0
//   d3 = a1*b2 + a2*b1
//   d4 = a2*b2
// Naming v0=a0*b0, v1=a1*b1, v2=a2*b2, each cross-term pair collapses to one
// extra mul the same way QuadExt::mul folds its own cross terms:
//   d1 = (a0+a1)*(b0+b1) - v0 - v1
//   d3 = (a1+a2)*(b1+b2) - v1 - v2
//   d2 = (a0+a2)*(b0+b2) - v0 - v2 + v1
// so d0..d4 cost 6 Base muls total (v0, v1, v2, and the three cross
// products) instead of 9. Reducing mod v^3 - XI folds v^3 -> XI and
// v^4 -> XI*v:
//   c0 = d0 + XI*d3
//   c1 = d1 + XI*d4 = d1 + XI*v2
//   c2 = d2
// The two XI multiplies are against the fixed constant XI, not one of the
// two operands being multiplied together, so (as with QuadExt::mul's BETA
// multiply) they're kept separate from the 6-mul count above.
impl<C: CubicExtConfig> Mul for CubicExt<C> {
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

        CubicExt::new(c0, c1, c2)
    }
}

// Lets CubicExt<C> itself serve as another extension's Base -- e.g. a
// further extension built on top of an Fp6-shaped tower.
impl<C: CubicExtConfig> Field for CubicExt<C> {
    fn zero() -> Self {
        CubicExt::new(C::Base::zero(), C::Base::zero(), C::Base::zero())
    }

    fn one() -> Self {
        CubicExt::new(C::Base::one(), C::Base::zero(), C::Base::zero())
    }

    fn inverse(self) -> Self {
        CubicExt::inverse(self)
    }
}

// Frobenius on Base[v]/(v^3-XI): x^p for x = c0 + c1*v + c2*v^2 expands,
// using that Frobenius is additive and multiplicative (char p), to
// c0.frobenius() + c1.frobenius()*v^p + c2.frobenius()*v^(2p). v^p and
// v^(2p) reduce to FROBENIUS_COEFF_C1*v and FROBENIUS_COEFF_C2*v^2 by
// CubicExtFrobeniusConfig's definitions, so this is c0.frobenius() +
// (FROBENIUS_COEFF_C1*c1.frobenius())*v + (FROBENIUS_COEFF_C2*
// c2.frobenius())*v^2 -- two Base multiplies, no exponentiation at
// runtime. Generic over any Base: Frobenius, same as QuadExt's Frobenius
// impl.
impl<C: CubicExtFrobeniusConfig> Frobenius for CubicExt<C>
where
    C::Base: Frobenius,
{
    fn frobenius(self) -> Self {
        CubicExt::new(
            self.c0.frobenius(),
            C::FROBENIUS_COEFF_C1 * self.c1.frobenius(),
            C::FROBENIUS_COEFF_C2 * self.c2.frobenius(),
        )
    }
}
