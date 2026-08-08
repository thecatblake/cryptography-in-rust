use std::marker::PhantomData;
use std::ops::{Add, Div, Mul, Neg, ShrAssign, Sub};

// The value should only satisfy this.
pub trait FpRepr: Copy + PartialEq + ShrAssign<usize> {
    const ZERO: Self;

    fn bit(&self, i: usize) -> bool;
}

impl FpRepr for u64 {
    const ZERO: Self = 0;

    fn bit(&self, i: usize) -> bool {
        (self >> i) & 1 == 1
    }
}

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
