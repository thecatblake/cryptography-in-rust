use field_core::Field;
use std::ops::{Add, Div, Mul, Rem, Sub};

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

    // Index of the highest-degree nonzero coefficient. Unlike coeffs.len(),
    // this ignores any trailing zero padding left behind by Add/Sub/Mul, so
    // it reflects the polynomial's true degree.
    fn degree(&self) -> Option<usize> {
        degree(&self.coeffs)
    }

    // Long division: repeatedly compare leading terms, divide self's
    // leading coefficient by g's, record that in the quotient, multiply g
    // by it and subtract from the remainder -- until the remainder's degree
    // drops below g's.
    pub fn div_rem(&self, g: &Self) -> (Self, Self) {
        let g_deg = g.degree().expect("division by the zero polynomial");
        let g_lead_inv = g.coeffs[g_deg].inverse();

        let mut remainder = self.coeffs.clone();
        let mut quotient = vec![F::zero(); remainder.len()];

        while let Some(r_deg) = degree(&remainder) {
            if r_deg < g_deg {
                break;
            }

            let shift = r_deg - g_deg;
            let term = remainder[r_deg] * g_lead_inv;
            quotient[shift] = term;

            for (i, &gc) in g.coeffs[..=g_deg].iter().enumerate() {
                remainder[shift + i] = remainder[shift + i] - term * gc;
            }
        }

        (DensePolynomial::new(quotient), DensePolynomial::new(remainder))
    }
}

// Free function (rather than a method) so it can be reused on plain slices
// without borrowing a DensePolynomial, e.g. while remainder is a bare Vec<F>
// mid-division.
fn degree<F: Field>(coeffs: &[F]) -> Option<usize> {
    coeffs.iter().rposition(|c| *c != F::zero())
}

impl<F: Field> Polynomial<F> for DensePolynomial<F> {
    // Horner's method: O(n) field muls instead of recomputing each x^i from scratch.
    fn evaluate(&self, x: F) -> F {
        self.coeffs.iter().rev().fold(F::zero(), |acc, &c| acc * x + c)
    }
}

fn zip_coeffs<'a, F: Field>(a: &'a [F], b: &'a [F]) -> impl Iterator<Item = (F, F)> + 'a {
    let len = a.len().max(b.len());

    (0..len).map(move |i| {
        let a = a.get(i).copied().unwrap_or(F::zero());
        let b = b.get(i).copied().unwrap_or(F::zero());

        (a, b)
    })
}

impl<F: Field> Add for DensePolynomial<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let coeffs = zip_coeffs(&self.coeffs, &rhs.coeffs).map(|(a, b)| a + b).collect();

        DensePolynomial::new(coeffs)
    }
}

impl<F: Field> Sub for DensePolynomial<F> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let coeffs = zip_coeffs(&self.coeffs, &rhs.coeffs).map(|(a, b)| a - b).collect();

        DensePolynomial::new(coeffs)
    }
}

impl<F: Field> Mul for DensePolynomial<F> {
    type Output = Self;

    // convolution: O(n*m) field muls over the coefficient lists.
    fn mul(self, rhs: Self) -> Self::Output {
        if self.coeffs.is_empty() || rhs.coeffs.is_empty() {
            return DensePolynomial::new(vec![]);
        }

        let mut coeffs = vec![F::zero(); self.coeffs.len() + rhs.coeffs.len() - 1];

        for (i, &a) in self.coeffs.iter().enumerate() {
            for (j, &b) in rhs.coeffs.iter().enumerate() {
                coeffs[i + j] = coeffs[i + j] + a * b;
            }
        }

        DensePolynomial::new(coeffs)
    }
}

impl<F: Field> Div for DensePolynomial<F> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self.div_rem(&rhs).0
    }
}

impl<F: Field> Rem for DensePolynomial<F> {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        self.div_rem(&rhs).1
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

    #[test]
    fn add_same_degree() {
        // (3 + 5x + 2x^2) + (1 + 2x + 9x^2) = 4 + 7x + 11x^2
        let p = DensePolynomial::new(vec![fe(3), fe(5), fe(2)]);
        let q = DensePolynomial::new(vec![fe(1), fe(2), fe(9)]);
        let sum = p + q;

        assert_eq!(sum.coeffs.len(), 3);
        assert_eq!(sum.coeffs[0].value, fe(4).value);
        assert_eq!(sum.coeffs[1].value, fe(7).value);
        assert_eq!(sum.coeffs[2].value, fe(11).value);
    }

    #[test]
    fn add_different_degree_pads_shorter_with_zero() {
        // (3 + 5x) + (1 + 2x + 9x^2) = 4 + 7x + 9x^2
        let p = DensePolynomial::new(vec![fe(3), fe(5)]);
        let q = DensePolynomial::new(vec![fe(1), fe(2), fe(9)]);
        let sum = p + q;

        assert_eq!(sum.coeffs.len(), 3);
        assert_eq!(sum.coeffs[0].value, fe(4).value);
        assert_eq!(sum.coeffs[1].value, fe(7).value);
        assert_eq!(sum.coeffs[2].value, fe(9).value);
    }

    #[test]
    fn sub_same_degree() {
        // (3 + 5x + 2x^2) - (1 + 2x + 9x^2) = 2 + 3x - 7x^2 (mod 17)
        let p = DensePolynomial::new(vec![fe(3), fe(5), fe(2)]);
        let q = DensePolynomial::new(vec![fe(1), fe(2), fe(9)]);
        let diff = p - q;

        assert_eq!(diff.coeffs.len(), 3);
        assert_eq!(diff.coeffs[0].value, fe(2).value);
        assert_eq!(diff.coeffs[1].value, fe(3).value);
        assert_eq!(diff.coeffs[2].value, fe(10).value);
    }

    #[test]
    fn sub_self_is_all_zero_coeffs() {
        let p = DensePolynomial::new(vec![fe(3), fe(5), fe(2)]);
        let coeffs_for_eval = p.coeffs.clone();
        let diff = DensePolynomial::new(coeffs_for_eval) - p;

        assert!(diff.coeffs.iter().all(|c| c.value == F17::zero().value));
    }

    proptest! {
        #[test]
        fn add_evaluate_is_homomorphic(
            p_coeffs in proptest::collection::vec(0u32..17, 0..8),
            q_coeffs in proptest::collection::vec(0u32..17, 0..8),
            x in 0u32..17,
        ) {
            let p_coeffs: Vec<F17> = p_coeffs.into_iter().map(fe).collect();
            let q_coeffs: Vec<F17> = q_coeffs.into_iter().map(fe).collect();
            let x = fe(x);

            let p = DensePolynomial::new(p_coeffs.clone());
            let q = DensePolynomial::new(q_coeffs.clone());
            let expected = DensePolynomial::new(p_coeffs).evaluate(x) + DensePolynomial::new(q_coeffs).evaluate(x);

            prop_assert_eq!((p + q).evaluate(x).value, expected.value);
        }

        #[test]
        fn sub_evaluate_is_homomorphic(
            p_coeffs in proptest::collection::vec(0u32..17, 0..8),
            q_coeffs in proptest::collection::vec(0u32..17, 0..8),
            x in 0u32..17,
        ) {
            let p_coeffs: Vec<F17> = p_coeffs.into_iter().map(fe).collect();
            let q_coeffs: Vec<F17> = q_coeffs.into_iter().map(fe).collect();
            let x = fe(x);

            let p = DensePolynomial::new(p_coeffs.clone());
            let q = DensePolynomial::new(q_coeffs.clone());
            let expected = DensePolynomial::new(p_coeffs).evaluate(x) - DensePolynomial::new(q_coeffs).evaluate(x);

            prop_assert_eq!((p - q).evaluate(x).value, expected.value);
        }
    }

    #[test]
    fn mul_by_zero_polynomial_is_empty() {
        let p = DensePolynomial::new(vec![fe(3), fe(5), fe(2)]);
        let zero = DensePolynomial::new(vec![]);
        let product = p * zero;

        assert!(product.coeffs.is_empty());
    }

    #[test]
    fn mul_example() {
        // (1 + 2x) * (3 + 4x) = 3 + 10x + 8x^2
        let p = DensePolynomial::new(vec![fe(1), fe(2)]);
        let q = DensePolynomial::new(vec![fe(3), fe(4)]);
        let product = p * q;

        assert_eq!(product.coeffs.len(), 3);
        assert_eq!(product.coeffs[0].value, fe(3).value);
        assert_eq!(product.coeffs[1].value, fe(10).value);
        assert_eq!(product.coeffs[2].value, fe(8).value);
    }

    #[test]
    fn mul_by_constant_scales_coeffs() {
        // (2) * (3 + 5x + 2x^2) = 6 + 10x + 4x^2
        let p = DensePolynomial::new(vec![fe(2)]);
        let q = DensePolynomial::new(vec![fe(3), fe(5), fe(2)]);
        let product = p * q;

        assert_eq!(product.coeffs.len(), 3);
        assert_eq!(product.coeffs[0].value, fe(6).value);
        assert_eq!(product.coeffs[1].value, fe(10).value);
        assert_eq!(product.coeffs[2].value, fe(4).value);
    }

    proptest! {
        #[test]
        fn mul_evaluate_is_homomorphic(
            p_coeffs in proptest::collection::vec(0u32..17, 0..6),
            q_coeffs in proptest::collection::vec(0u32..17, 0..6),
            x in 0u32..17,
        ) {
            let p_coeffs: Vec<F17> = p_coeffs.into_iter().map(fe).collect();
            let q_coeffs: Vec<F17> = q_coeffs.into_iter().map(fe).collect();
            let x = fe(x);

            let p = DensePolynomial::new(p_coeffs.clone());
            let q = DensePolynomial::new(q_coeffs.clone());
            let expected = DensePolynomial::new(p_coeffs).evaluate(x) * DensePolynomial::new(q_coeffs).evaluate(x);

            prop_assert_eq!((p * q).evaluate(x).value, expected.value);
        }

        #[test]
        fn mul_degree_matches_sum_of_degrees(
            p_coeffs in proptest::collection::vec(1u32..17, 1..6),
            q_coeffs in proptest::collection::vec(1u32..17, 1..6),
        ) {
            let p_coeffs: Vec<F17> = p_coeffs.into_iter().map(fe).collect();
            let q_coeffs: Vec<F17> = q_coeffs.into_iter().map(fe).collect();
            let expected_len = p_coeffs.len() + q_coeffs.len() - 1;

            let p = DensePolynomial::new(p_coeffs);
            let q = DensePolynomial::new(q_coeffs);

            prop_assert_eq!((p * q).coeffs.len(), expected_len);
        }
    }

    fn trimmed_values(coeffs: &[F17]) -> Vec<u32> {
        let mut values: Vec<u32> = coeffs.iter().map(|c| c.value).collect();
        while values.last() == Some(&0) {
            values.pop();
        }
        values
    }

    #[test]
    fn div_rem_exact_division() {
        // x^3 - 1 = (x - 1)(x^2 + x + 1), so dividing by (x - 1) leaves no remainder.
        let p = DensePolynomial::new(vec![fe(16), fe(0), fe(0), fe(1)]); // -1 + x^3
        let g = DensePolynomial::new(vec![fe(16), fe(1)]); // -1 + x
        let (q, r) = p.div_rem(&g);

        assert_eq!(trimmed_values(&q.coeffs), vec![1, 1, 1]);
        assert!(trimmed_values(&r.coeffs).is_empty());
    }

    #[test]
    fn div_rem_non_monic_divisor() {
        // 2x^2 + 4x = (2x + 4) * x, remainder 0.
        let p = DensePolynomial::new(vec![fe(0), fe(4), fe(2)]);
        let g = DensePolynomial::new(vec![fe(4), fe(2)]);
        let (q, r) = p.div_rem(&g);

        assert_eq!(trimmed_values(&q.coeffs), vec![0, 1]);
        assert!(trimmed_values(&r.coeffs).is_empty());
    }

    #[test]
    fn div_rem_with_nonzero_remainder() {
        // (x^2 + 1) / x = x remainder 1.
        let p = DensePolynomial::new(vec![fe(1), fe(0), fe(1)]);
        let g = DensePolynomial::new(vec![fe(0), fe(1)]);
        let (q, r) = p.div_rem(&g);

        assert_eq!(trimmed_values(&q.coeffs), vec![0, 1]);
        assert_eq!(trimmed_values(&r.coeffs), vec![1]);
    }

    #[test]
    fn div_rem_dividend_degree_less_than_divisor() {
        let p = DensePolynomial::new(vec![fe(3), fe(5)]);
        let g = DensePolynomial::new(vec![fe(1), fe(2), fe(9)]);
        let (q, r) = p.div_rem(&g);

        assert!(trimmed_values(&q.coeffs).is_empty());
        assert_eq!(trimmed_values(&r.coeffs), trimmed_values(&p.coeffs));
    }

    #[test]
    #[should_panic(expected = "division by the zero polynomial")]
    fn div_rem_by_zero_polynomial_panics() {
        let p = DensePolynomial::new(vec![fe(1), fe(2)]);
        let zero: DensePolynomial<F17> = DensePolynomial::new(vec![]);

        p.div_rem(&zero);
    }

    proptest! {
        #[test]
        fn div_rem_reconstructs_dividend(
            p_coeffs in proptest::collection::vec(0u32..17, 0..8),
            g_coeffs in proptest::collection::vec(1u32..17, 1..5),
        ) {
            let p_coeffs: Vec<F17> = p_coeffs.into_iter().map(fe).collect();
            let g_coeffs: Vec<F17> = g_coeffs.into_iter().map(fe).collect();

            let p = DensePolynomial::new(p_coeffs.clone());
            let g = DensePolynomial::new(g_coeffs.clone());
            let (q, r) = p.div_rem(&g);

            let g_deg = degree(&g_coeffs).unwrap();
            prop_assert!(r.degree().is_none_or(|d| d < g_deg));

            let reconstructed = DensePolynomial::new(q.coeffs.clone()) * DensePolynomial::new(g_coeffs) + DensePolynomial::new(r.coeffs.clone());
            prop_assert_eq!(trimmed_values(&reconstructed.coeffs), trimmed_values(&p_coeffs));
        }

        #[test]
        fn div_operator_matches_div_rem_quotient(
            p_coeffs in proptest::collection::vec(0u32..17, 0..8),
            g_coeffs in proptest::collection::vec(1u32..17, 1..5),
        ) {
            let p_coeffs: Vec<F17> = p_coeffs.into_iter().map(fe).collect();
            let g_coeffs: Vec<F17> = g_coeffs.into_iter().map(fe).collect();

            let p = DensePolynomial::new(p_coeffs.clone());
            let g = DensePolynomial::new(g_coeffs.clone());
            let expected = DensePolynomial::new(p_coeffs.clone()).div_rem(&DensePolynomial::new(g_coeffs.clone())).0;

            prop_assert_eq!(trimmed_values(&(p / g).coeffs), trimmed_values(&expected.coeffs));
        }

        #[test]
        fn rem_operator_matches_div_rem_remainder(
            p_coeffs in proptest::collection::vec(0u32..17, 0..8),
            g_coeffs in proptest::collection::vec(1u32..17, 1..5),
        ) {
            let p_coeffs: Vec<F17> = p_coeffs.into_iter().map(fe).collect();
            let g_coeffs: Vec<F17> = g_coeffs.into_iter().map(fe).collect();

            let p = DensePolynomial::new(p_coeffs.clone());
            let g = DensePolynomial::new(g_coeffs.clone());
            let expected = DensePolynomial::new(p_coeffs.clone()).div_rem(&DensePolynomial::new(g_coeffs.clone())).1;

            prop_assert_eq!(trimmed_values(&(p % g).coeffs), trimmed_values(&expected.coeffs));
        }
    }
}
