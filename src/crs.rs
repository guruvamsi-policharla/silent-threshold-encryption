use crate::utils::lagrange_poly;
use ark_ec::{pairing::Pairing, ScalarMul};
use ark_ec::{PrimeGroup, VariableBaseMSM};
use ark_ff::{Field, PrimeField};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_std::rand::Rng;
use ark_std::{One, UniformRand, Zero};

#[derive(Clone, Debug)]
pub struct CRS<E: Pairing> {
    pub n: usize,
    pub powers_of_g: Vec<E::G1Affine>,
    pub powers_of_h: Vec<E::G2Affine>,

    // preprocessed lagrange polynomials
    pub li: Vec<E::G1>,
    pub li_minus0: Vec<E::G1>,
    pub li_x: Vec<E::G1>,
    pub li_lj_z: Vec<Vec<E::G1>>,

    // preprocessed lagrange polynomials in g2 (only needed for verifying hints)
    pub li_g2: Vec<E::G2>,
    pub li_minus0_g2: Vec<E::G2>,
    pub li_x_g2: Vec<E::G2>,
    pub li_lj_z_g2: Vec<Vec<E::G2>>,
}

impl<E: Pairing> CRS<E> {
    pub fn new(n: usize, rng: &mut impl Rng) -> Self {
        // powers of tau
        let tau = E::ScalarField::rand(rng);

        let mut powers_of_tau = vec![E::ScalarField::one()];

        let mut cur = tau;
        for _ in 0..=n {
            powers_of_tau.push(cur);
            cur *= &tau;
        }

        let powers_of_g = E::G1::generator().batch_mul(&powers_of_tau[0..n + 1]);
        let powers_of_h = E::G2::generator().batch_mul(&powers_of_tau[0..n + 1]);

        // lagrange powers
        let mut li_evals: Vec<E::ScalarField> = vec![E::ScalarField::zero(); n];
        let mut li_evals_minus0: Vec<E::ScalarField> = vec![E::ScalarField::zero(); n];
        let mut li_evals_x: Vec<E::ScalarField> = vec![E::ScalarField::zero(); n];

        let tau2_inv: <E as Pairing>::ScalarField = (tau * tau).inverse().unwrap();
        for i in 0..n {
            let li = lagrange_poly(n, i);
            li_evals[i] = li.evaluate(&tau);

            li_evals_minus0[i] = (li_evals[i] - li.coeffs[0]) * tau;

            li_evals_x[i] = li_evals_minus0[i] * tau2_inv;
        }

        let z_eval = tau.pow(&[n as u64]) - E::ScalarField::one();
        let z_eval_inv = z_eval.inverse().unwrap();

        let mut li = vec![E::G1::zero(); n];
        let mut li_g2 = vec![E::G2::zero(); n];

        for i in 0..n {
            li[i] = E::G1::generator() * li_evals[i];
            li_g2[i] = E::G2::generator() * li_evals[i];
        }

        let mut li_minus0 = vec![E::G1::zero(); n];
        let mut li_minus0_g2 = vec![E::G2::zero(); n];

        for i in 0..n {
            li_minus0[i] = E::G1::generator() * li_evals_minus0[i];
            li_minus0_g2[i] = E::G2::generator() * li_evals_minus0[i];
        }

        let mut li_x = vec![E::G1::zero(); n];
        let mut li_x_g2 = vec![E::G2::zero(); n];

        for i in 0..n {
            li_x[i] = E::G1::generator() * li_evals_x[i];
            li_x_g2[i] = E::G2::generator() * li_evals_x[i];
        }

        let mut li_lj_z = vec![vec![E::G1::zero(); n]; n];
        let mut li_lj_z_g2 = vec![vec![E::G2::zero(); n]; n];

        for i in 0..n {
            for j in 0..n {
                li_lj_z[i][j] = if i == j {
                    E::G1::generator() * ((li_evals[i] * li_evals[i] - li_evals[i]) * z_eval_inv)
                } else {
                    E::G1::generator() * (li_evals[i] * li_evals[j] * z_eval_inv)
                };

                li_lj_z_g2[i][j] = if i == j {
                    E::G2::generator() * ((li_evals[i] * li_evals[i] - li_evals[i]) * z_eval_inv)
                } else {
                    E::G2::generator() * (li_evals[i] * li_evals[j] * z_eval_inv)
                };
            }
        }

        Self {
            n,
            powers_of_g,
            powers_of_h,

            li,
            li_minus0,
            li_x,
            li_lj_z,

            li_g2,
            li_minus0_g2,
            li_x_g2,
            li_lj_z_g2,
        }
    }

    pub fn commit_g1(&self, coeffs: &Vec<E::ScalarField>) -> E::G1 {
        assert!(
            coeffs.len() <= self.powers_of_g.len(),
            "Too many coefficients for the given powers of tau"
        );

        let plain_coeffs = coeffs.iter().map(|c| c.into_bigint()).collect::<Vec<_>>();
        <E::G1 as VariableBaseMSM>::msm_bigint(
            &self.powers_of_g[..coeffs.len()],
            plain_coeffs.as_slice(),
        )
    }

    pub fn commit_g2(&self, coeffs: &Vec<E::ScalarField>) -> E::G2 {
        assert!(
            coeffs.len() <= self.powers_of_g.len(),
            "Too many coefficients for the given powers of tau"
        );

        let plain_coeffs = coeffs.iter().map(|c| c.into_bigint()).collect::<Vec<_>>();
        <E::G2 as VariableBaseMSM>::msm_bigint(
            &self.powers_of_h[..coeffs.len()],
            plain_coeffs.as_slice(),
        )
    }

    pub fn compute_opening_proof(
        &self,
        coeffs: &Vec<E::ScalarField>,
        point: &E::ScalarField,
    ) -> E::G1 {
        let polynomial = DensePolynomial::from_coefficients_slice(&coeffs);
        let eval = polynomial.evaluate(point);

        let mut numerator = polynomial.clone();
        numerator.coeffs[0] -= eval;

        let divisor = DensePolynomial::from_coefficients_vec(vec![
            E::ScalarField::zero() - point,
            E::ScalarField::one(),
        ]);
        let witness_polynomial = &numerator / &divisor;

        self.commit_g1(&witness_polynomial.coeffs)
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Bls12_381 as E;
    use ark_bls12_381::Fr as F;
    use ark_bls12_381::G1Projective as G1;
    use ark_bls12_381::G2Projective as G2;
    use ark_ec::pairing::Pairing;
    use ark_ec::PrimeGroup;
    use ark_poly::EvaluationDomain;
    use ark_poly::Polynomial;
    use ark_poly::Radix2EvaluationDomain;
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
    use ark_std::UniformRand;
    use ark_std::Zero;

    #[test]
    fn test_sumcheck() {
        // A(X).B(X) = \sum_i A(i).B(i) + X * Q_x(X) + Z(X) * Q_Z(X)
        let rng = &mut ark_std::test_rng();

        let n = 1 << 5;
        let domain = Radix2EvaluationDomain::<F>::new(n).unwrap();

        // sample n random evals
        let a_evals = (0..n).map(|_| F::rand(rng)).collect::<Vec<_>>();
        let b_evals = (0..n).map(|_| F::rand(rng)).collect::<Vec<_>>();
        let mut s = F::zero();
        for i in 0..n {
            s += a_evals[i] * b_evals[i];
        }

        let a_coeffs = domain.ifft(&a_evals);
        let b_coeffs = domain.ifft(&b_evals);

        let a_poly = DensePolynomial::from_coefficients_vec(a_coeffs);
        let b_poly = DensePolynomial::from_coefficients_vec(b_coeffs);

        let c_poly = &a_poly * &b_poly;

        println!("a_poly deg: {}", a_poly.degree());
        println!("b_poly deg: {}", b_poly.degree());
        println!("c_poly deg: {}", c_poly.degree());

        let (qz, rem) = c_poly.divide_by_vanishing_poly(domain);
        println!("qz deg: {}", qz.degree());
        println!("rem deg: {}", rem.degree());

        assert_eq!(s / F::from(n as u64), rem.evaluate(&F::zero()));
    }

    #[test]
    fn test_kzg() {
        // A(X).B(X) = \sum_i A(i).B(i) + X * Q_x(X) + Z(X) * Q_Z(X)
        let rng = &mut ark_std::test_rng();

        let n = 1 << 3;
        let crs = crate::crs::CRS::<E>::new(n, rng);

        // sample n random coeffs
        let coeffs = (0..n).map(|_| F::rand(rng)).collect::<Vec<_>>();
        let com = crs.commit_g1(&coeffs);

        let point = F::rand(rng);
        let eval = DensePolynomial::from_coefficients_slice(&coeffs).evaluate(&point);

        let pi = crs.compute_opening_proof(&coeffs, &point);

        let lhs = E::pairing(com + (G1::generator() * (-eval)), G2::generator());
        let rhs = E::pairing(pi, crs.powers_of_h[1] - (G2::generator() * point));
        assert_eq!(lhs, rhs);
    }
}
