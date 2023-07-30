use ark_ff::{FftField, Field};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, Evaluations, Polynomial, Radix2EvaluationDomain,
};

//returns t(X) = X^n - 1
// pub fn compute_vanishing_poly<F: Field>(n: usize) -> DensePolynomial<F> {
//     let mut coeffs = vec![];
//     for i in 0..n+1 {
//         if i == 0 {
//             coeffs.push(F::zero() - F::one()); // -1
//         } else if i == n {
//             coeffs.push(F::one()); // X^n
//         } else {
//             coeffs.push(F::zero());
//         }
//     }
//     DensePolynomial { coeffs }
// }

// 1 at omega^i and 0 elsewhere on domain {omega^i}_{i \in [n]}
pub fn lagrange_poly<F: FftField>(n: usize, i: usize) -> DensePolynomial<F> {
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

// returns t(X) = X
pub fn compute_x_monomial<F: Field>() -> DensePolynomial<F> {
    let mut coeffs = vec![];
    coeffs.push(F::zero()); // 0
    coeffs.push(F::one()); // X
    DensePolynomial { coeffs }
}

// returns t(X) = c
pub fn compute_constant_poly<F: Field>(c: &F) -> DensePolynomial<F> {
    let mut coeffs = vec![];
    coeffs.push(c.clone()); // c
    DensePolynomial { coeffs }
}

//computes f(ωx)
pub fn poly_domain_mult_ω<F: Field>(f: &DensePolynomial<F>, ω: &F) -> DensePolynomial<F> {
    let mut new_poly = f.clone();
    for i in 1..(f.degree() + 1) {
        //we don't touch the zeroth coefficient
        let ω_pow_i: F = ω.pow([i as u64]);
        new_poly.coeffs[i] = new_poly.coeffs[i] * ω_pow_i;
    }
    new_poly
}

//computes c . f(x), for some constnt c
pub fn poly_eval_mult_c<F: Field>(f: &DensePolynomial<F>, c: &F) -> DensePolynomial<F> {
    let mut new_poly = f.clone();
    for i in 0..(f.degree() + 1) {
        new_poly.coeffs[i] = new_poly.coeffs[i] * c.clone();
    }
    new_poly
}
