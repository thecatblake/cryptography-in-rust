use field_core::Field;

pub trait Polynomial<F: Field> {
    fn evaluate(&self, x: F) -> F;
}

pub struct DensePolynomial<F: Field> {
    // coeffs[i] is the coefficient of x^i, so coeffs[0] is the constant term.
    pub coeffs: Vec<F>,
}

impl<F: Field> DensePolynomial<F> {
    pub fn new(coeffs: Vec<F>) -> Self {
        DensePolynomial { coeffs }
    }
}

impl<F: Field> Polynomial<F> for DensePolynomial<F> {
    // Horner's method: O(n) field muls instead of recomputing each x^i from scratch.
    fn evaluate(&self, x: F) -> F {
        self.coeffs.iter().rev().fold(F::zero(), |acc, &c| acc * x + c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use field_core::{Fp, WideArithmeticBackend, WideFieldConfig};
    use proptest::prelude::*;

    struct Mod17;
    impl WideFieldConfig for Mod17 {
        type Repr = u32;

        const MODULUS: u32 = 17;

        fn mul(a: u32, b: u32) -> u32 {
            ((a as u64 * b as u64) % 17) as u32
        }
    }

    type B17 = WideArithmeticBackend<Mod17>;
    type F17 = Fp<B17>;

    fn fe(v: u32) -> F17 {
        F17::new(v % 17)
    }

    fn eval_naive(coeffs: &[F17], x: F17) -> F17 {
        let mut result = F17::zero();
        let mut power = F17::one();

        for &c in coeffs {
            result = result + c * power;
            power = power * x;
        }

        result
    }

    #[test]
    fn evaluate_constant_polynomial() {
        let p = DensePolynomial::new(vec![fe(5)]);
        assert_eq!(p.evaluate(fe(3)).value, fe(5).value);
    }

    #[test]
    fn evaluate_empty_polynomial_is_zero() {
        let p: DensePolynomial<F17> = DensePolynomial::new(vec![]);
        assert_eq!(p.evaluate(fe(9)).value, F17::zero().value);
    }

    #[test]
    fn evaluate_at_zero_is_constant_term() {
        let p = DensePolynomial::new(vec![fe(4), fe(2), fe(9)]);
        assert_eq!(p.evaluate(F17::zero()).value, fe(4).value);
    }

    #[test]
    fn evaluate_matches_naive_sum() {
        // p(x) = 3 + 5x + 2x^2
        let p = DensePolynomial::new(vec![fe(3), fe(5), fe(2)]);
        let x = fe(6);

        assert_eq!(p.evaluate(x).value, eval_naive(&p.coeffs, x).value);
    }

    proptest! {
        #[test]
        fn evaluate_matches_naive_sum_proptest(
            coeffs in proptest::collection::vec(0u32..17, 0..8),
            x in 0u32..17,
        ) {
            let coeffs: Vec<F17> = coeffs.into_iter().map(fe).collect();
            let p = DensePolynomial::new(coeffs.clone());
            let x = fe(x);

            prop_assert_eq!(p.evaluate(x).value, eval_naive(&coeffs, x).value);
        }
    }
}
