use ark_ff::{FftField, Field};
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Evaluations, Polynomial,
    Radix2EvaluationDomain,
};

// 1 at omega^i and 0 elsewhere on domain {omega^i}_{i \in [n]}
pub fn lagrange_poly<F: FftField>(n: usize, i: usize) -> DensePolynomial<F> {
    debug_assert!(i < n);
    //todo: check n is a power of 2
    let mut evals = vec![];
    for j in 0..n {
        let l_of_x: u64 = if i == j { 1 } else { 0 };
        evals.push(F::from(l_of_x));
    }

    //powers of nth root of unity
    let domain = Radix2EvaluationDomain::<F>::new(n).unwrap();
    let eval_form = Evaluations::from_vec_and_domain(evals, domain);
    //interpolated polynomial over the n points
    eval_form.interpolate()
}

/// interpolates a polynomial which is zero on points and 1 at the point 0
/// todo: use faster interpolation
pub fn interp_mostly_zero<F: Field>(points: &Vec<F>) -> DensePolynomial<F> {
    if points.is_empty() {
        // threshold=n
        return DensePolynomial::from_coefficients_vec(vec![F::one()]);
    }

    let mut interp = DensePolynomial::from_coefficients_vec(vec![F::one()]);
    for &point in points {
        interp = interp.naive_mul(&DensePolynomial::from_coefficients_vec(vec![
            -point,
            F::one(),
        ]));
    }

    let scale = interp.evaluate(&F::zero());
    interp = &interp * (F::one() / scale);

    interp
}

// /// interpolates a polynomial where evaluations on points are zero and the polynomial evaluates to 1 at the point 1
// pub fn compute_vanishing_poly(points: &Vec<ScalarField>) -> DensePolynomial {
//     let mut monomials = Vec::new();
//     for i in 0..points.len() {
//         monomials.push(DensePolynomial::from_coeffs(
//             HostSlice::from_slice(&vec![ScalarField::zero() - points[i], ScalarField::one()]),
//             2,
//         ));
//     }

//     // assert that points.len() is a power of 2
//     assert_eq!(
//         points.len().count_ones(),
//         1,
//         "Implementation demands that n-t is a power of 2. Currently: {}",
//         points.len()
//     );

//     let mut chunk_size = points.len() / 2;
//     while chunk_size > 0 {
//         for i in 0..chunk_size {
//             monomials[i] = &monomials[i] * &monomials[i + chunk_size];
//         }
//         chunk_size = chunk_size / 2;
//     }

//     let scale = monomials[0].eval(&ScalarField::one());
//     let res = &monomials[0] * &scale.inv();

//     res
// }
