use std::marker::PhantomData;
use std::ops::{Add, Div, Mul, Neg, Shr, ShrAssign, Sub};

mod quad_ext;
pub use quad_ext::{QuadExt, QuadExtConfig};

mod cubic_ext;
pub use cubic_ext::{CubicExt, CubicExtConfig};

mod fp2;
pub use fp2::{Fp2, Fp2Config};

mod fp6;
pub use fp6::{Fp6, Fp6Config};

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

pub trait FpBackend {
    type Repr: FpRepr;

    const MODULUS: Self::Repr;

    fn add(a: Self::Repr, b: Self::Repr) -> Self::Repr;
    fn sub(a: Self::Repr, b: Self::Repr) -> Self::Repr;
    fn mul(a: Self::Repr, b: Self::Repr) -> Self::Repr;
    fn neg(a: Self::Repr) -> Self::Repr;
    fn inverse(a: Self::Repr) -> Self::Repr;
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
// on it aren't tied to primitive integers only.
pub trait WideInt: FpRepr {
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

    // Fermat's little theorem: a^(p-2) == a^-1 (mod p).
    fn inverse(a: T::Repr) -> T::Repr {
        assert!(a != T::Repr::ZERO, "cannot invert zero in a field");

        let mut result = T::Repr::from_u8(1);
        let mut base = a;
        let mut e = T::Repr::narrow(Self::MODULUS.widen() - T::Repr::from_u8(2).widen());

        while e != T::Repr::ZERO {
            if e.bit(0) {
                result = Self::mul(result, base);
            }
            base = Self::square(base);
            e >>= 1;
        }

        result
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

    // Fermat's little theorem, carried out entirely in Montgomery form:
    // Montgomery multiplication is a ring isomorphism, so repeated
    // Self::mul/Self::square on Montgomery-form values computes
    // (a^e)_mont directly, with no to_mont/from_mont round-trip needed.
    fn inverse(a: T::Repr) -> T::Repr {
        assert!(a != T::Repr::ZERO, "cannot invert zero in a field");

        let mut result = Self::one();
        let mut base = a;
        let mut e = T::Repr::narrow(Self::MODULUS.widen() - T::Repr::from_u8(2).widen());

        while e != T::Repr::ZERO {
            if e.bit(0) {
                result = Self::mul(result, base);
            }
            base = Self::square(base);
            e >>= 1;
        }

        result
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
pub trait Field: Copy + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Neg<Output = Self> {
    // Defaults to self*self; Fp<B> overrides it to go through
    // FpBackend::square, which some backends implement more cheaply than a
    // general multiply.
    fn square(self) -> Self {
        self * self
    }

    fn inverse(self) -> Self;
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
    fn square(self) -> Self {
        Fp::square(self)
    }

    fn inverse(self) -> Self {
        Fp::inverse(self)
    }
}
