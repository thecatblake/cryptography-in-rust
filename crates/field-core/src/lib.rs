use std::marker::PhantomData;
use std::ops::{Add, Div, Mul, Neg, Shr, ShrAssign, Sub};

mod quad_ext;
pub use quad_ext::{QuadExt, QuadExtConfig, QuadExtFrobeniusConfig};

mod cubic_ext;
pub use cubic_ext::{CubicExt, CubicExtConfig, CubicExtFrobeniusConfig};

mod fp2;
pub use fp2::{Fp2, Fp2Config};

mod fp6;
pub use fp6::{Fp6, Fp6Config};

mod fp12;
pub use fp12::{Fp12, Fp12Config};

// The value should only satisfy this.
pub trait FpRepr: Copy + PartialEq + ShrAssign<usize> {
    const ZERO: Self;

    fn bit(&self, i: usize) -> bool;
}

macro_rules! impl_fp_repr {
    ($repr:ty) => {
        impl FpRepr for $repr {
            const ZERO: Self = 0;

            fn bit(&self, i: usize) -> bool {
                (self >> i) & 1 == 1
            }
        }
    };
}

impl_fp_repr!(u8);
impl_fp_repr!(u16);
impl_fp_repr!(u32);
impl_fp_repr!(u64);
impl_fp_repr!(u128);

// Repr types that support extended-GCD-based modular inversion. Sub/mul
// wrap on overflow (two's-complement semantics), matching bigint::Uint<N>'s
// Add/Sub impls -- the extended Euclidean algorithm's Bezout coefficient
// goes "negative" mid-algorithm and relies on that wraparound, recovered
// via is_negative's top-bit check (see gcd_inverse below).
pub trait EuclideanRepr: FpRepr {
    const ONE: Self;

    fn wrapping_sub(self, rhs: Self) -> Self;
    fn wrapping_mul(self, rhs: Self) -> Self;
    fn div_rem(self, rhs: Self) -> (Self, Self);
    fn is_negative(self) -> bool;
}

macro_rules! impl_euclidean_repr {
    ($repr:ty) => {
        impl EuclideanRepr for $repr {
            const ONE: Self = 1;

            fn wrapping_sub(self, rhs: Self) -> Self {
                <$repr>::wrapping_sub(self, rhs)
            }

            fn wrapping_mul(self, rhs: Self) -> Self {
                <$repr>::wrapping_mul(self, rhs)
            }

            fn div_rem(self, rhs: Self) -> (Self, Self) {
                (self / rhs, self % rhs)
            }

            fn is_negative(self) -> bool {
                self.bit(<$repr>::BITS as usize - 1)
            }
        }
    };
}

impl_euclidean_repr!(u8);
impl_euclidean_repr!(u16);
impl_euclidean_repr!(u32);
impl_euclidean_repr!(u64);
impl_euclidean_repr!(u128);

// Extended-Euclidean-algorithm modular inverse: computes a^-1 mod modulus.
// Only tracks the Bezout coefficient for `a`, not `modulus`'s -- that's all
// a modular inverse needs, unlike a full extended-gcd. Sub/mul wrap
// (two's-complement semantics), so the coefficient is free to go
// "negative" mid-algorithm; is_negative's top-bit check plus a final
// negate-and-subtract-from-modulus recovers the canonical positive residue.
pub fn gcd_inverse<R: EuclideanRepr>(a: R, modulus: R) -> R {
    let mut old_r = a;
    let mut r = modulus;
    let mut old_s = R::ONE;
    let mut s = R::ZERO;

    while r != R::ZERO {
        let (q, rem) = old_r.div_rem(r);

        old_r = r;
        r = rem;

        let new_s = old_s.wrapping_sub(q.wrapping_mul(s));
        old_s = s;
        s = new_s;
    }

    if old_s.is_negative() {
        modulus.wrapping_sub(R::ZERO.wrapping_sub(old_s))
    } else {
        old_s
    }
}

pub trait FpBackend {
    // EuclideanRepr (not just FpRepr) so every backend's Repr can drive the
    // default `inverse` below -- see gcd_inverse's doc comment for why
    // that's the default. A hypothetical Repr that couldn't support it
    // would still need to override inverse, but would also need its own
    // div_rem/wrapping_sub/wrapping_mul to satisfy this bound in the first
    // place, so in practice every integer-like Repr qualifies for free.
    type Repr: FpRepr + EuclideanRepr;

    const MODULUS: Self::Repr;

    fn add(a: Self::Repr, b: Self::Repr) -> Self::Repr;
    fn sub(a: Self::Repr, b: Self::Repr) -> Self::Repr;
    fn mul(a: Self::Repr, b: Self::Repr) -> Self::Repr;
    fn neg(a: Self::Repr) -> Self::Repr;
    // Extended-GCD-based inversion by default (see gcd_inverse) -- far
    // cheaper than Fermat's little theorem for wide moduli (e.g. ~450x for
    // a 256-bit field, since GCD's cost tracks the operand's bit-pattern
    // rather than a full modular exponentiation). Backends with a cheaper
    // option (e.g. WideArithmeticBackend/MontWideBackend's fermat_inverse
    // for machine-word fields) override this.
    fn inverse(a: Self::Repr) -> Self::Repr {
        assert!(a != Self::Repr::ZERO, "cannot invert zero in a field");
        gcd_inverse(a, Self::MODULUS)
    }
    // The multiplicative identity in this backend's representation
    // (plain 1 for DefaultBackend, R mod MODULUS for MontBackend).
    fn one() -> Self::Repr;

    // Defaults to mul(a, a); backends override this when a dedicated
    // squaring routine is cheaper than a general multiply.
    fn square(a: Self::Repr) -> Self::Repr {
        Self::mul(a, a)
    }
}

// An integer with a strictly wider scratch type to add/multiply into
// without overflow before reducing back down. Not "native machine int"
// specific by design -- a bigint type with a double-width counterpart
// (e.g. U256 with U512) can satisfy this same contract, so backends built
// on it aren't tied to primitive integers only. Also EuclideanRepr so
// WideArithmeticBackend/MontWideBackend's Repr satisfies FpBackend::Repr's
// bound even though both backends override `inverse` with fermat_inverse
// instead of using the EuclideanRepr-driven default.
pub trait WideInt: FpRepr + EuclideanRepr {
    type Wide: Copy + PartialOrd + Add<Output = Self::Wide> + Sub<Output = Self::Wide> + Shr<usize, Output = Self::Wide>;

    // Bit width of Self (so R = 2^BITS for Montgomery-style backends).
    const BITS: usize;

    fn widen(self) -> Self::Wide;
    fn narrow(wide: Self::Wide) -> Self;
    fn wide_mul(self, other: Self) -> Self::Wide;
    fn from_u8(v: u8) -> Self;
}

// Only pairs with a strictly wider native integer type can implement this
// (there's no built-in "next size up" past u128), so each pair is spelled
// out once here rather than derived generically.
macro_rules! impl_wide_int {
    ($repr:ty => $wide:ty) => {
        impl WideInt for $repr {
            type Wide = $wide;

            const BITS: usize = <$repr>::BITS as usize;

            fn widen(self) -> $wide {
                self as $wide
            }

            fn narrow(wide: $wide) -> $repr {
                wide as $repr
            }

            fn wide_mul(self, other: $repr) -> $wide {
                self as $wide * other as $wide
            }

            fn from_u8(v: u8) -> $repr {
                v as $repr
            }
        }
    };
}

impl_wide_int!(u8 => u16);
impl_wide_int!(u16 => u32);
impl_wide_int!(u32 => u64);
impl_wide_int!(u64 => u128);

// a^(MODULUS-2) mod p via square-and-multiply -- Fermat's little theorem.
// Same loop shape as Fp::pow below, but operating on the raw Repr since
// FpBackend::inverse is called before an Fp<B> value exists to invoke
// .pow() on. Opt-in alternative to the EuclideanRepr-driven GCD default:
// cheap and simple for a machine-word modulus (small fixed iteration
// count), unlike the ~450x-slower Fermat exponentiation a wide (e.g.
// 256-bit) modulus would need -- see gcd_inverse's doc comment.
pub fn fermat_inverse<B: FpBackend>(a: B::Repr) -> B::Repr
where
    B::Repr: WideInt,
{
    let mut result = B::one();
    let mut base = a;
    let mut e = B::Repr::narrow(B::MODULUS.widen() - B::Repr::from_u8(2).widen());

    while e != B::Repr::ZERO {
        if e.bit(0) {
            result = B::mul(result, base);
        }
        base = B::square(base);
        e >>= 1;
    }

    result
}

// Small fields (representable in a single machine word) share the same
// add/sub/neg/inverse shape; only the multiplication reduction differs
// enough per-field to be worth hand-optimizing (e.g. Goldilocks' epsilon
// trick). Implementors only need to supply the representation width, the
// modulus, and that one routine.
pub trait WideFieldConfig {
    type Repr: WideInt;

    const MODULUS: Self::Repr;

    fn mul(a: Self::Repr, b: Self::Repr) -> Self::Repr;
}

pub struct WideArithmeticBackend<T: WideFieldConfig>(PhantomData<T>);

impl<T: WideFieldConfig> FpBackend for WideArithmeticBackend<T> {
    type Repr = T::Repr;

    const MODULUS: T::Repr = T::MODULUS;

    fn add(a: T::Repr, b: T::Repr) -> T::Repr {
        let sum = a.widen() + b.widen();
        let modulus = Self::MODULUS.widen();

        T::Repr::narrow(if sum >= modulus { sum - modulus } else { sum })
    }

    fn sub(a: T::Repr, b: T::Repr) -> T::Repr {
        let a = a.widen();
        let b = b.widen();

        T::Repr::narrow(if a >= b { a - b } else { a + Self::MODULUS.widen() - b })
    }

    fn mul(a: T::Repr, b: T::Repr) -> T::Repr {
        T::mul(a, b)
    }

    fn neg(a: T::Repr) -> T::Repr {
        if a == T::Repr::ZERO { T::Repr::ZERO } else { T::Repr::narrow(Self::MODULUS.widen() - a.widen()) }
    }

    // Fermat's little theorem instead of the EuclideanRepr-driven GCD
    // default -- cheap and simple at machine-word size (see
    // fermat_inverse's doc comment).
    fn inverse(a: T::Repr) -> T::Repr {
        assert!(a != T::Repr::ZERO, "cannot invert zero in a field");
        fermat_inverse::<Self>(a)
    }

    fn one() -> T::Repr {
        T::Repr::from_u8(1)
    }
}

// A field represented via Montgomery (REDC) reduction instead of a native
// reduction trick: values are stored as a*R mod MODULUS ("Montgomery
// form"), and REDC undoes that scaling as a side effect of reducing a
// product mod MODULUS, using only shifts, adds, and multiplies -- no
// division.
//
// Only safe when MODULUS < R/2 (R = 2^Repr::BITS), i.e. MODULUS uses at
// most BITS-1 bits: REDC's intermediate `t + m*MODULUS` is bounded by
// 2*R*MODULUS, which must fit inside Wide (2*BITS bits) without an extra
// guard bit. A modulus using the full width of Repr (e.g. a 256-bit prime
// stored in a 256-bit Repr) needs one more bit of scratch than Wide
// provides, which this backend doesn't allocate.
pub trait MontFieldConfig {
    type Repr: WideInt;

    const MODULUS: Self::Repr;
    const R2: Self::Repr;
    const N_PRIME: Self::Repr;
}

pub struct MontWideBackend<T: MontFieldConfig>(PhantomData<T>);

impl<T: MontFieldConfig> MontWideBackend<T> {
    fn redc(t: <T::Repr as WideInt>::Wide) -> T::Repr {
        let t_low = T::Repr::narrow(t);
        let m = T::Repr::narrow(t_low.wide_mul(T::N_PRIME));

        // t + m*MODULUS is guaranteed divisible by R by construction of m.
        let sum = t + m.wide_mul(T::MODULUS);
        let quotient = sum >> T::Repr::BITS;

        let modulus = T::MODULUS.widen();

        T::Repr::narrow(if quotient >= modulus { quotient - modulus } else { quotient })
    }

    pub fn to_mont(a: T::Repr) -> T::Repr {
        Self::redc(a.wide_mul(T::R2))
    }

    pub fn from_mont(a: T::Repr) -> T::Repr {
        Self::redc(a.widen())
    }
}

// Const-evaluable counterpart to MontWideBackend::to_mont, one per
// primitive width (mirrors impl_fp_repr!/impl_wide_int!'s per-width
// codegen above). A fully generic `to_mont<T: MontFieldConfig>` can't be a
// const fn on stable Rust: internally it calls WideInt's widen/wide_mul/
// narrow through generic trait dispatch, and const trait dispatch isn't
// stable yet. Written directly against a concrete primitive pair instead
// -- same REDC steps as MontWideBackend::to_mont/redc, just with `as`
// casts standing in for widen/narrow/wide_mul -- so it *is* const, just
// not generic. Lets a MontFieldConfig/Fp2Config implementor write BETA (or
// any other Montgomery-form constant) starting from its plain/canonical
// value, with the conversion happening at compile time instead of via a
// runtime to_mont() call or an offline-computed magic number.
macro_rules! impl_const_to_mont {
    ($name:ident, $repr:ty, $wide:ty) => {
        pub const fn $name(value: $repr, r2: $repr, n_prime: $repr, modulus: $repr) -> $repr {
            let t = (value as $wide) * (r2 as $wide);
            let t_low = t as $repr;
            let m = ((t_low as $wide) * (n_prime as $wide)) as $repr;

            // t + m*modulus is guaranteed divisible by R by construction of m.
            let sum = t + (m as $wide) * (modulus as $wide);
            let quotient = sum >> <$repr>::BITS;

            let modulus_wide = modulus as $wide;

            (if quotient >= modulus_wide { quotient - modulus_wide } else { quotient }) as $repr
        }
    };
}

impl_const_to_mont!(to_mont_u8, u8, u16);
impl_const_to_mont!(to_mont_u16, u16, u32);
impl_const_to_mont!(to_mont_u32, u32, u64);
impl_const_to_mont!(to_mont_u64, u64, u128);

// Const-evaluable modular exponentiation on Montgomery-form values, one per
// primitive width (mirrors impl_const_to_mont! above, for the same reason:
// a fully generic const version can't exist on stable Rust, since it would
// need const trait dispatch through WideInt). Given `base` already in
// Montgomery form, computes base^exp via square-and-multiply, using a local
// mont_mul (same REDC steps as impl_const_to_mont!'s inline reduction) so
// the loop can reduce after both the squaring and the conditional multiply
// without duplicating those steps. Lets a Fp2Config implementor write
// FROBENIUS_COEFF as pow_mont_uNN(<Self as Fp2Config>::BETA.value,
// (Self::MODULUS - 1) / 2, Self::R2, Self::N_PRIME, Self::MODULUS) --
// BETA^((p-1)/2), entirely at compile time.
macro_rules! impl_const_pow_mont {
    ($name:ident, $repr:ty, $wide:ty) => {
        pub const fn $name(base: $repr, exp: $repr, r2: $repr, n_prime: $repr, modulus: $repr) -> $repr {
            const fn mont_mul(a: $repr, b: $repr, n_prime: $repr, modulus: $repr) -> $repr {
                let t = (a as $wide) * (b as $wide);
                let t_low = t as $repr;
                let m = ((t_low as $wide) * (n_prime as $wide)) as $repr;

                // t + m*modulus is guaranteed divisible by R by construction of m.
                let sum = t + (m as $wide) * (modulus as $wide);
                let quotient = sum >> <$repr>::BITS;

                let modulus_wide = modulus as $wide;

                (if quotient >= modulus_wide { quotient - modulus_wide } else { quotient }) as $repr
            }

            // Montgomery-form 1, i.e. REDC(1*R2) = REDC(R2) -- same value
            // MontWideBackend::one produces at runtime -- seeds the
            // square-and-multiply accumulator.
            let mut result = mont_mul(1, r2, n_prime, modulus);
            let mut b = base;
            let mut e = exp;

            while e != 0 {
                if e & 1 == 1 {
                    result = mont_mul(result, b, n_prime, modulus);
                }
                b = mont_mul(b, b, n_prime, modulus);
                e >>= 1;
            }

            result
        }
    };
}

impl_const_pow_mont!(pow_mont_u8, u8, u16);
impl_const_pow_mont!(pow_mont_u16, u16, u32);
impl_const_pow_mont!(pow_mont_u32, u32, u64);
impl_const_pow_mont!(pow_mont_u64, u64, u128);

// Const-evaluable modular exponentiation of a Montgomery-form Fp2 element
// (a0 + a1*u, u^2 = beta), one per primitive width -- the Fp2-level
// counterpart to impl_const_pow_mont! above. Needed because an Fp6Config's
// Frobenius coefficient (Fp2Config::FROBENIUS_COEFF is a base-field scalar,
// but Fp6Config::FROBENIUS_COEFF_C1/C2 are Fp2 elements, XI^((p-1)/3) and
// XI^(2(p-1)/3)) can't be computed by pow_mont_uNN alone: it needs
// quadratic-extension multiplication, not just base-field multiplication,
// and QuadExt::mul goes through Field's regular (non-const) trait methods.
// So this hand-rolls the same Karatsuba shape as QuadExt::mul (qmul below)
// plus mod-p add/sub, entirely in terms of the primitive-width REDC steps
// already used by impl_const_pow_mont!, and runs square-and-multiply over
// that. Lets a Fp6Config implementor write FROBENIUS_COEFF_C1 as
// pow_mont_fp2_uNN(Self::XI.c0.value, Self::XI.c1.value,
// (Self::MODULUS - 1) / 3, <Self as Fp2Config>::BETA.value, Self::R2,
// Self::N_PRIME, Self::MODULUS) (see field_core::pow_mont_fp2_u32 et al.),
// wrapping the returned (c0, c1) pair in Fp2 { c0: Fp::new(c0), c1:
// Fp::new(c1) } -- FROBENIUS_COEFF_C2 uses twice that exponent.
macro_rules! impl_const_pow_mont_fp2 {
    ($name:ident, $repr:ty, $wide:ty) => {
        pub const fn $name(
            a0: $repr,
            a1: $repr,
            exp: $repr,
            beta: $repr,
            r2: $repr,
            n_prime: $repr,
            modulus: $repr,
        ) -> ($repr, $repr) {
            const fn mont_mul(a: $repr, b: $repr, n_prime: $repr, modulus: $repr) -> $repr {
                let t = (a as $wide) * (b as $wide);
                let t_low = t as $repr;
                let m = ((t_low as $wide) * (n_prime as $wide)) as $repr;

                // t + m*modulus is guaranteed divisible by R by construction of m.
                let sum = t + (m as $wide) * (modulus as $wide);
                let quotient = sum >> <$repr>::BITS;

                let modulus_wide = modulus as $wide;

                (if quotient >= modulus_wide { quotient - modulus_wide } else { quotient }) as $repr
            }

            const fn mont_add(a: $repr, b: $repr, modulus: $repr) -> $repr {
                let sum = (a as $wide) + (b as $wide);
                let modulus_wide = modulus as $wide;

                (if sum >= modulus_wide { sum - modulus_wide } else { sum }) as $repr
            }

            const fn mont_sub(a: $repr, b: $repr, modulus: $repr) -> $repr {
                let a = a as $wide;
                let b = b as $wide;
                let modulus_wide = modulus as $wide;

                (if a >= b { a - b } else { a + modulus_wide - b }) as $repr
            }

            // Quadratic-extension Montgomery multiplication: same Karatsuba
            // shape as QuadExt::mul (c0 = v0 + beta*v1, c1 = (a0+a1)(b0+b1)
            // - v0 - v1), built from mont_mul/mont_add/mont_sub above.
            const fn qmul(
                a0: $repr,
                a1: $repr,
                b0: $repr,
                b1: $repr,
                beta: $repr,
                n_prime: $repr,
                modulus: $repr,
            ) -> ($repr, $repr) {
                let v0 = mont_mul(a0, b0, n_prime, modulus);
                let v1 = mont_mul(a1, b1, n_prime, modulus);

                let c0 = mont_add(v0, mont_mul(beta, v1, n_prime, modulus), modulus);
                let cross =
                    mont_mul(mont_add(a0, a1, modulus), mont_add(b0, b1, modulus), n_prime, modulus);
                let c1 = mont_sub(mont_sub(cross, v0, modulus), v1, modulus);

                (c0, c1)
            }

            // Montgomery-form (1, 0), the Fp2 multiplicative identity --
            // seeds the square-and-multiply accumulator.
            let mut result0 = mont_mul(1, r2, n_prime, modulus);
            let mut result1: $repr = 0;

            let mut base0 = a0;
            let mut base1 = a1;
            let mut e = exp;

            while e != 0 {
                if e & 1 == 1 {
                    let (c0, c1) = qmul(result0, result1, base0, base1, beta, n_prime, modulus);
                    result0 = c0;
                    result1 = c1;
                }
                let (s0, s1) = qmul(base0, base1, base0, base1, beta, n_prime, modulus);
                base0 = s0;
                base1 = s1;
                e >>= 1;
            }

            (result0, result1)
        }
    };
}

impl_const_pow_mont_fp2!(pow_mont_fp2_u8, u8, u16);
impl_const_pow_mont_fp2!(pow_mont_fp2_u16, u16, u32);
impl_const_pow_mont_fp2!(pow_mont_fp2_u32, u32, u64);
impl_const_pow_mont_fp2!(pow_mont_fp2_u64, u64, u128);

impl<T: MontFieldConfig> FpBackend for MontWideBackend<T> {
    type Repr = T::Repr;

    const MODULUS: T::Repr = T::MODULUS;

    fn add(a: T::Repr, b: T::Repr) -> T::Repr {
        let sum = a.widen() + b.widen();
        let modulus = Self::MODULUS.widen();

        T::Repr::narrow(if sum >= modulus { sum - modulus } else { sum })
    }

    fn sub(a: T::Repr, b: T::Repr) -> T::Repr {
        let a = a.widen();
        let b = b.widen();

        T::Repr::narrow(if a >= b { a - b } else { a + Self::MODULUS.widen() - b })
    }

    fn mul(a: T::Repr, b: T::Repr) -> T::Repr {
        Self::redc(a.wide_mul(b))
    }

    fn neg(a: T::Repr) -> T::Repr {
        if a == T::Repr::ZERO { T::Repr::ZERO } else { T::Repr::narrow(Self::MODULUS.widen() - a.widen()) }
    }

    // Fermat's little theorem instead of the EuclideanRepr-driven GCD
    // default, same as WideArithmeticBackend (see fermat_inverse's doc
    // comment). Carried out entirely in Montgomery form for free: Montgomery
    // multiplication is a ring isomorphism, so fermat_inverse's repeated
    // Self::mul/Self::square on Montgomery-form values computes (a^e)_mont
    // directly, with no to_mont/from_mont round-trip needed.
    fn inverse(a: T::Repr) -> T::Repr {
        assert!(a != T::Repr::ZERO, "cannot invert zero in a field");
        fermat_inverse::<Self>(a)
    }

    fn one() -> T::Repr {
        Self::to_mont(T::Repr::from_u8(1))
    }
}

// Field abstracts over "a ring QuadExt/CubicExt can be built on top of": the
// four ring operations plus squaring and inversion. Fp<B> is the base case
// below; QuadExt<C> and CubicExt<C> each implement it too (in terms of
// their own Base: Field), which is what lets extensions stack on top of
// each other -- e.g. a cubic extension of a quadratic one -- without
// QuadExt/CubicExt needing to know anything Fp-specific about Base.
//
// PartialEq is a supertrait, not an add-on bound callers reach for: a
// field's elements have plain, total equality (no NaN-like case to carve
// out, same reasoning FpRepr's own PartialEq supertrait above rests on),
// and every implementor here (Fp<B>, QuadExt<C>, CubicExt<C>) already
// provides it. Requiring it here means downstream generic code (e.g.
// elliptic_curve::AffinePoint) gets equality for free instead of threading
// a `where C::Field: PartialEq` bound through every impl that needs it.
pub trait Field:
    Copy + PartialEq + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Neg<Output = Self>
{
    // The additive and multiplicative identities. Not consts: Fp<B>'s ONE is
    // FpBackend::one(), which for MontBackend is R mod MODULUS rather than a
    // literal 1, and that conversion isn't available as a const fn over a
    // generic B. Needed by callers that must produce a field element from
    // nothing -- e.g. twisted_edwards::AffinePoint's identity (0,1) -- since
    // unlike short_weierstrass's point-at-infinity flag, that identity is a
    // genuine affine point with real coordinates.
    fn zero() -> Self;
    fn one() -> Self;

    // Defaults to self*self; Fp<B> overrides it to go through
    // FpBackend::square, which some backends implement more cheaply than a
    // general multiply.
    fn square(self) -> Self {
        self * self
    }

    fn inverse(self) -> Self;
}

// The Frobenius endomorphism x -> x^p, where p is the field's own prime
// characteristic. On Fp itself this is the identity (Fermat's little
// theorem: a^p == a for every a in Fp) -- the interesting cases are
// extension fields, where x^p permutes basis coefficients instead of
// fixing everything. Kept as its own trait rather than a defaulted method
// on Field: an "identity by default" default would be silently wrong for
// any extension that forgot to override it, and unlike square/inverse
// there's no formula generic over an arbitrary Field -- each extension
// needs its own extension-specific coefficient (see Fp2Config::
// FROBENIUS_COEFF).
pub trait Frobenius: Field {
    fn frobenius(self) -> Self;
}

pub struct Fp<B: FpBackend> {
    pub value: B::Repr,
    _marker: PhantomData<B>,
}

// Derived impls would require B: Clone/Copy, but only B::Repr needs to be
// (FpRepr's Copy supertrait already guarantees that), so implement by hand.
impl<B: FpBackend> Clone for Fp<B> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<B: FpBackend> Copy for Fp<B> {}

// Both backends keep every value in a single canonical representation
// (DefaultBackend: reduced mod p; MontBackend: Montgomery form, whose
// to_mont/from_mont map is a bijection), so comparing the stored
// representations directly is equivalent to comparing the field elements
// they denote -- no reduction or from_mont round-trip needed first.
impl<B: FpBackend> PartialEq for Fp<B> {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl<B: FpBackend> Fp<B> {
    pub const fn new(value: B::Repr) -> Self {
        Fp { value, _marker: PhantomData }
    }

    pub fn inverse(self) -> Self {
        Fp::new(B::inverse(self.value))
    }

    pub fn square(self) -> Self {
        Fp::new(B::square(self.value))
    }

    // Square-and-multiply modular exponentiation: O(log exp) field muls
    // instead of O(exp).
    pub fn pow(self, exp: B::Repr) -> Self {
        let mut result = B::one();
        let mut base = self.value;
        let mut e = exp;

        while e != B::Repr::ZERO {
            if e.bit(0) {
                result = B::mul(result, base);
            }
            base = B::square(base);
            e >>= 1;
        }

        Fp::new(result)
    }
}

impl<B: FpBackend> Add for Fp<B> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Fp::new(B::add(self.value, rhs.value))
    }
}

impl<B: FpBackend> Sub for Fp<B> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Fp::new(B::sub(self.value, rhs.value))
    }
}

impl<B: FpBackend> Mul for Fp<B> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Fp::new(B::mul(self.value, rhs.value))
    }
}

impl<B: FpBackend> Div for Fp<B> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.inverse()
    }
}

impl<B: FpBackend> Neg for Fp<B> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Fp::new(B::neg(self.value))
    }
}

impl<B: FpBackend> Field for Fp<B> {
    fn zero() -> Self {
        Fp::new(B::Repr::ZERO)
    }

    fn one() -> Self {
        Fp::new(B::one())
    }

    fn square(self) -> Self {
        Fp::square(self)
    }

    fn inverse(self) -> Self {
        Fp::inverse(self)
    }
}

impl<B: FpBackend> Frobenius for Fp<B> {
    // Fermat's little theorem: a^p == a (mod p) for every a in Fp.
    fn frobenius(self) -> Self {
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    // Ported from bigint::math's retired extended_gcd/mod_inv tests (now
    // that gcd_inverse lives here and is generic over EuclideanRepr instead
    // of hardcoded to Uint<N>) -- same known-answer vectors, run over u32.
    #[test]
    fn mod_inv_3_11() {
        assert_eq!(gcd_inverse(3u32, 11), 4);
    }

    #[test]
    fn mod_inv_one() {
        assert_eq!(gcd_inverse(1u32, 7), 1);
    }

    #[test]
    fn mod_inv_rsa_example() {
        // e = 17, phi = 3120 -> d = 2753 (17 * 2753 = 46801 = 15*3120 + 1)
        assert_eq!(gcd_inverse(17u32, 3120), 2753);
    }

    #[test]
    fn mod_inv_reduces_a_greater_than_n() {
        // 40 mod 7 = 5, and 5 * 3 = 15 = 1 mod 7
        assert_eq!(gcd_inverse(40u32, 7), 3);
    }

    fn gcd_u32(mut a: u32, mut b: u32) -> u32 {
        while b != 0 {
            let t = b;
            b = a % b;
            a = t;
        }
        a
    }

    proptest! {
        #[test]
        fn gcd_inverse_is_multiplicative_inverse(a in 1u32.., n in 2u32..) {
            prop_assume!(gcd_u32(a, n) == 1);

            let inv = gcd_inverse(a % n, n);

            prop_assert!(inv < n);
            prop_assert_eq!((a as u64 * inv as u64) % (n as u64), 1);
        }
    }

    // Toy field (mod 97, u32) purely to compare the two inversion
    // algorithms against each other -- field-core can't depend on
    // babybear/goldilocks/mersenne31 (they depend on it, not the other way
    // around), so this stands in for "a real WideArithmeticBackend".
    struct Mod97;
    impl WideFieldConfig for Mod97 {
        type Repr = u32;

        const MODULUS: u32 = 97;

        fn mul(a: u32, b: u32) -> u32 {
            ((a as u64 * b as u64) % 97) as u32
        }
    }
    type B97 = WideArithmeticBackend<Mod97>;

    proptest! {
        // WideArithmeticBackend::inverse (Fermat, via fermat_inverse) must
        // agree with the EuclideanRepr-driven GCD default it opts out of.
        #[test]
        fn fermat_inverse_matches_gcd_inverse(a in 1u32..97) {
            prop_assert_eq!(B97::inverse(a), gcd_inverse(a, B97::MODULUS));
        }
    }
}
