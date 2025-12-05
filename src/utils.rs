use crate::polynomial::{DensePolynomial, Evaluations, Radix2EvaluationDomain};
use blstrs::Scalar;
use ff::{BatchInvert, Field};

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

/// Compute all Lagrange basis polynomials over an n-point radix-2 domain.
/// Returns polynomials L_i such that L_i(omega^j) = 1 if i == j else 0.
pub fn lagrange_polys(n: usize) -> Vec<DensePolynomial> {
    assert!(n.is_power_of_two(), "domain size must be power of two");
    let domain = Radix2EvaluationDomain::new(n).expect("valid radix-2 domain");
    let omega = domain.group_gen;
    let omega_inv = omega.invert().unwrap();
    let n_scalar = Scalar::from(n as u64);

    // precompute omega^{-i} for i in 0..n
    let mut omega_inv_pows = Vec::with_capacity(n);
    let mut cur = Scalar::ONE;
    for _ in 0..n {
        omega_inv_pows.push(cur);
        cur *= omega_inv;
    }

    // compute 1 / (n * omega^{-i}) for all i via batch inversion
    let mut denom_invs: Vec<Scalar> = omega_inv_pows.iter().map(|w| *w * n_scalar).collect();
    denom_invs.iter_mut().batch_invert();

    let mut result = Vec::with_capacity(n);
    for (omega_i_inv, denom_inv) in omega_inv_pows.iter().zip(denom_invs.iter()) {
        // coefficients: x^{n-1} + omega_i x^{n-2} + ... + omega_i^{n-1}
        // stored low-to-high degree
        let mut coeffs = Vec::with_capacity(n);
        let mut power = *omega_i_inv; // omega_i^{-1}
        for _ in 0..n {
            coeffs.push(power * denom_inv);
            power *= *omega_i_inv; // multiply by omega_i^{-1} each step
        }

        result.push(DensePolynomial::from_coefficients_vec(coeffs));
    }

    result
}

/// interpolates a polynomial when all evaluations except at points[0] are zero
/// todo: check that multiplication is fast as one polynomial is shorter
pub fn interp_mostly_zero(eval: Scalar, points: &[Scalar]) -> DensePolynomial {
    if points.is_empty() {
        // threshold=n
        return DensePolynomial::from_coefficients_vec(vec![Scalar::ONE]);
    }

    // Build product \prod_{i>0}(x - points[i]) iteratively on coefficients
    let mut coeffs = vec![Scalar::ONE];
    for &point in points.iter().skip(1) {
        let neg_point = -point;
        coeffs.push(Scalar::ZERO);
        for i in (0..coeffs.len() - 1).rev() {
            let (left, right) = coeffs.split_at_mut(i + 1);
            let coef = &mut left[i];
            let next = &mut right[0];
            *next += *coef;
            *coef *= neg_point;
        }
    }

    // Evaluate at points[0] (Horner)
    let mut scale = *coeffs.last().unwrap();
    for c in coeffs.iter().rev().skip(1) {
        scale = scale * points[0] + c;
    }

    let scale_inv = scale.invert().unwrap();
    for c in coeffs.iter_mut() {
        *c *= eval * scale_inv;
    }

    DensePolynomial::from_coefficients_vec(coeffs)
}
