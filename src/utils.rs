use ark_ff::{FftField, Field};
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Evaluations, Polynomial,
    Radix2EvaluationDomain,
};

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

// interpolate a polynomial given evaluations on a set of points
pub fn lagrange_interpolate<F: Field>(evals: &Vec<F>, points: &Vec<F>) -> DensePolynomial<F> {
    // compute lagrange polynomials over the domain "points"
    let mut lagrange_polys = vec![];
    for i in 0..points.len() {
        let mut num = DensePolynomial::from_coefficients_vec(vec![F::one()]);
        let mut den = F::one();
        for j in 0..points.len() {
            if i == j {
                continue;
            }
            let x_xj = DensePolynomial::from_coefficients_vec(vec![-points[j], F::one()]);
            num = num.naive_mul(&x_xj);
            // let numerator = DensePolynomial::from_coefficients_vec(vec![-points[j], F::one()]);
            den *= points[i] - points[j];
            // lagrage_poly = num/den
        }
        den.inverse_in_place();
        lagrange_polys.push(&num * den);
    }

    // compute the interpolation
    let mut interp = DensePolynomial::from_coefficients_vec(vec![F::zero()]);
    for i in 0..points.len() {
        interp = interp + &lagrange_polys[i] * evals[i];
    }

    interp
}

// interpolates a polynomial when all evaluations except at points[0] are zero
// todo: check that multiplication is fast as one polynomial is shorter
pub fn interp_mostly_zero<F: Field>(eval: F, points: &Vec<F>) -> DensePolynomial<F> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ec::pairing::Pairing;
    use ark_std::{One, Zero};

    type E = ark_bls12_381::Bls12_381;
    type F = <E as Pairing>::ScalarField;

    #[test]
    fn interpolate_test() {
        let evals = vec![F::one(); 5];
        let points = vec![
            F::from(1u32),
            F::from(2u32),
            F::from(3u32),
            F::from(4u32),
            F::from(5u32),
        ];

        let interp = lagrange_interpolate(&evals, &points);
        assert_eq!(interp.coeffs, vec![F::one()]);

        let evals = vec![F::one(), F::zero(), F::zero(), F::zero(), F::zero()];
        let points = vec![
            F::from(1u32),
            F::from(2u32),
            F::from(3u32),
            F::from(4u32),
            F::from(5u32),
        ];

        let interp = lagrange_interpolate(&evals, &points);
        let fast_interp = interp_mostly_zero(F::one(), &points);

        assert_eq!(interp, fast_interp);

        assert_eq!(interp.evaluate(&F::from(1u32)), F::one());
        assert_eq!(interp.evaluate(&F::from(2u32)), F::zero());
        assert_eq!(interp.evaluate(&F::from(3u32)), F::zero());
        assert_eq!(interp.evaluate(&F::from(4u32)), F::zero());
        assert_eq!(interp.evaluate(&F::from(5u32)), F::zero());
    }
}
