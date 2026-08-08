use std::marker::PhantomData;
use std::ops::{Add, Div, Mul, Neg, ShrAssign, Sub};

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

// A machine integer that's compact enough to be a native field element
// (u32, u64, ...), plus a wider integer type it can be widened into to do
// add/sub without overflow before reducing back down.
pub trait NativeInt: FpRepr {
    type Wide: Copy + PartialOrd + Add<Output = Self::Wide> + Sub<Output = Self::Wide>;

    fn widen(self) -> Self::Wide;
    fn narrow(wide: Self::Wide) -> Self;
    fn from_u8(v: u8) -> Self;
}

// Only pairs with a strictly wider native integer type can implement this
// (there's no built-in "next size up" past u128), so each pair is spelled
// out once here rather than derived generically.
macro_rules! impl_native_int {
    ($repr:ty => $wide:ty) => {
        impl NativeInt for $repr {
            type Wide = $wide;

            fn widen(self) -> $wide {
                self as $wide
            }

            fn narrow(wide: $wide) -> $repr {
                wide as $repr
            }

            fn from_u8(v: u8) -> $repr {
                v as $repr
            }
        }
    };
}

impl_native_int!(u8 => u16);
impl_native_int!(u16 => u32);
impl_native_int!(u32 => u64);
impl_native_int!(u64 => u128);

// Small native fields (representable in a single u32/u64) share the same
// add/sub/neg/inverse shape; only the multiplication reduction differs
// enough per-field to be worth hand-optimizing (e.g. Goldilocks' epsilon
// trick vs BabyBear's plain `%`). Implementors only need to supply the
// representation width, the modulus, and that one routine.
pub trait NativeFieldConfig {
    type Repr: NativeInt;

    const MODULUS: Self::Repr;

    fn mul(a: Self::Repr, b: Self::Repr) -> Self::Repr;
}

pub struct NativeArithmeticBackend<T: NativeFieldConfig>(PhantomData<T>);

impl<T: NativeFieldConfig> FpBackend for NativeArithmeticBackend<T> {
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
    pub fn new(value: B::Repr) -> Self {
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
