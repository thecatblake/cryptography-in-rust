use std::ops::{Add, Mul, Neg, Sub};

use crate::{Field, Frobenius};

// QuadExtConfig describes a quadratic extension Base[u] / (u^2 - BETA) for
// some base field Base. Fp2Config (Base = Fp<MontWideBackend<C>>) is the
// only instantiation used elsewhere in this crate, but nothing here is tied
// to Fp specifically: any Field works as Base, including QuadExt/CubicExt
// themselves, so towers of extensions (e.g. a quadratic extension of a
// cubic one) fall out for free.
pub trait QuadExtConfig: Sized {
    type Base: Field;

    // Must be a quadratic non-residue in Base, or u^2 - BETA factors and
    // the extension collapses to Base x Base instead of being a field.
    const BETA: Self::Base;
}

// QuadExtFrobeniusConfig supplies the one extra piece of data Frobenius
// needs beyond QuadExtConfig: BETA^((p-1)/2) in Base, so that u^p =
// u*(u^2)^((p-1)/2) reduces to FROBENIUS_COEFF*u (see the Frobenius impl
// below). Kept separate from QuadExtConfig itself, rather than folding
// FROBENIUS_COEFF into it directly, so that computing this coefficient --
// which needs the base field's actual characteristic, and in practice
// Montgomery-form compile-time exponentiation -- isn't forced on every
// QuadExtConfig implementor. Fp2Config is this crate's only implementor
// (fp2.rs blanket-impls this trait for any C: Fp2Config); Fp12's
// Fp12Marker doesn't implement it yet, so Fp12 doesn't implement Frobenius
// yet either -- that's a matter of adding the impl, not changing this
// trait or the Frobenius impl below.
pub trait QuadExtFrobeniusConfig: QuadExtConfig {
    const FROBENIUS_COEFF: Self::Base;
}

// QuadExt = Base[u] / (u^2 - BETA), elements represented as c0 + c1*u.
pub struct QuadExt<C: QuadExtConfig> {
    pub c0: C::Base,
    pub c1: C::Base,
}

// Derived impls would require C: Clone/Copy, but only C::Base needs to be
// (Field's Copy supertrait already guarantees that), so implement by hand
// -- same reasoning as Fp<B>'s hand-written Clone/Copy in lib.rs.
impl<C: QuadExtConfig> Clone for QuadExt<C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<C: QuadExtConfig> Copy for QuadExt<C> {}

// No `where C::Base: PartialEq` needed: C::Base: Field, and Field now
// carries PartialEq as a supertrait.
impl<C: QuadExtConfig> PartialEq for QuadExt<C> {
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1
    }
}

impl<C: QuadExtConfig> QuadExt<C> {
    pub fn new(c0: C::Base, c1: C::Base) -> Self {
        QuadExt { c0, c1 }
    }

    // N(a) = a * conjugate(a) = (c0+c1*u)(c0-c1*u) = c0^2 - c1^2*u^2, and
    // u^2 == BETA by definition of QuadExt, so this collapses to a
    // base-field element. `inverse` below is defined in terms of it.
    pub fn norm(self) -> C::Base {
        self.c0.square() - C::BETA * self.c1.square()
    }

    // Multiplicative inverse. For a = c0 + c1*u,
    // a^-1 = conjugate(a) / norm(a) = (c0 - c1*u) * norm(a)^-1. The
    // denominator is the same c0^2 - BETA*c1^2 as `norm`, computed in place
    // here rather than via self.norm() since it's needed before the
    // conjugate is built.
    pub fn inverse(self) -> Self {
        let denom_inv = (self.c0.square() - C::BETA * self.c1.square()).inverse();

        QuadExt::new(self.c0 * denom_inv, -(self.c1 * denom_inv))
    }
}

impl<C: QuadExtConfig> Add for QuadExt<C> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        QuadExt::new(self.c0 + rhs.c0, self.c1 + rhs.c1)
    }
}

impl<C: QuadExtConfig> Sub for QuadExt<C> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        QuadExt::new(self.c0 - rhs.c0, self.c1 - rhs.c1)
    }
}

impl<C: QuadExtConfig> Neg for QuadExt<C> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        QuadExt::new(-self.c0, -self.c1)
    }
}

// Karatsuba: schoolbook (a0+a1*u)(b0+b1*u) needs 4 base-field muls
// (a0*b0, a1*b1, a0*b1, a1*b0). Naming v0 = a0*b0 and v1 = a1*b1, the two
// cross terms' sum a0*b1 + a1*b0 equals (a0+a1)*(b0+b1) - v0 - v1, so only
// one more mul -- not two -- is needed to get both of them at once:
//   c0 = v0 + BETA*v1
//   c1 = (a0+a1)*(b0+b1) - v0 - v1
// 3 base-field muls total (v0, v1, and the cross-term product); the BETA
// multiply is separate since BETA is a fixed constant, not one of the two
// operands being multiplied together.
impl<C: QuadExtConfig> Mul for QuadExt<C> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let v0 = self.c0 * rhs.c0;
        let v1 = self.c1 * rhs.c1;

        let c0 = v0 + C::BETA * v1;
        let c1 = (self.c0 + self.c1) * (rhs.c0 + rhs.c1) - v0 - v1;

        QuadExt::new(c0, c1)
    }
}

// Lets QuadExt<C> itself serve as another extension's Base -- e.g.
// CubicExtConfig::Base = QuadExt<C> for a cubic-over-quadratic tower.
impl<C: QuadExtConfig> Field for QuadExt<C> {
    fn inverse(self) -> Self {
        QuadExt::inverse(self)
    }
}

// Frobenius on Base[u]/(u^2-BETA): x^p for x = c0 + c1*u expands, using
// that Frobenius is additive and multiplicative (both hold in any
// characteristic-p ring), to c0^p + c1^p*u^p = c0.frobenius() +
// c1.frobenius()*u^p. u^p reduces to FROBENIUS_COEFF*u by
// QuadExtFrobeniusConfig's definition, so this is c0.frobenius() +
// (FROBENIUS_COEFF*c1.frobenius())*u -- one Base multiply, no
// exponentiation at runtime. Generic over any Base: Frobenius (not just
// Fp), so this also covers towers built on top of QuadExt, once their
// QuadExtFrobeniusConfig is supplied.
impl<C: QuadExtFrobeniusConfig> Frobenius for QuadExt<C>
where
    C::Base: Frobenius,
{
    fn frobenius(self) -> Self {
        QuadExt::new(self.c0.frobenius(), C::FROBENIUS_COEFF * self.c1.frobenius())
    }
}
