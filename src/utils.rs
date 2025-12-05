use crate::polynomial::{DensePolynomial, Evaluations, Radix2EvaluationDomain};
use blstrs::Scalar;
use ff::Field;

// 1 at omega^i and 0 elsewhere on domain {omega^i}_{i \in [n]}
pub fn lagrange_poly(n: usize, i: usize) -> DensePolynomial {
    debug_assert!(i < n);
    //todo: check n is a power of 2
    let mut evals = vec![];
    for j in 0..n {
        let l_of_x = if i == j { Scalar::ONE } else { Scalar::ZERO };
        evals.push(l_of_x);
    }

    //powers of nth root of unity
    let domain = Radix2EvaluationDomain::new(n).unwrap();
    let eval_form = Evaluations::from_vec_and_domain(evals, domain);
    //interpolated polynomial over the n points
    eval_form.interpolate()
}

/// interpolates a polynomial when all evaluations except at points[0] are zero
/// todo: check that multiplication is fast as one polynomial is shorter
pub fn interp_mostly_zero(eval: Scalar, points: &Vec<Scalar>) -> DensePolynomial {
    if points.is_empty() {
        // threshold=n
        return DensePolynomial::from_coefficients_vec(vec![Scalar::ONE]);
    }

    let mut interp = DensePolynomial::from_coefficients_vec(vec![Scalar::ONE]);
    for &point in &points[1..] {
        interp = interp.naive_mul(&DensePolynomial::from_coefficients_vec(vec![
            -point,
            Scalar::ONE,
        ]));
    }

    let scale = interp.evaluate(&points[0]);
    interp = &interp * (eval * scale.invert().unwrap());

    interp
}
