use ark_ff::Field;
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};

/// interpolates a polynomial when all evaluations except at points[0] are zero
/// todo: check that multiplication is fast as one polynomial is shorter
pub fn interp_mostly_zero<F: Field>(eval: F, points: &Vec<F>) -> DensePolynomial<F> {
    if points.is_empty() {
        // threshold=n
        return DensePolynomial::from_coefficients_vec(vec![F::one()]);
    }

    let mut interp = DensePolynomial::from_coefficients_vec(vec![F::one()]);
    for &point in &points[1..] {
        interp = interp.naive_mul(&DensePolynomial::from_coefficients_vec(vec![
            -point,
            F::one(),
        ]));
    }

    let scale = interp.evaluate(&points[0]);
    interp = &interp * (eval / scale);

    interp
}
