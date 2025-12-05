//adapted from https://github.com/arkworks-rs/poly-commit/blob/master/src/kzg10/mod.rs

use blstrs::{G1Affine, G1Projective, G2Affine, G2Projective, Gt, Scalar};
use ff::Field;
use group::{Curve, Group};
use pairing::Engine;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

use crate::polynomial::DensePolynomial;

pub struct KZG10;

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct PowersOfTau {
    /// Group elements of the form `{ \tau^i G }`, where `i` ranges from 0 to `degree`.
    pub powers_of_g: Vec<G1Affine>,
    /// Group elements of the form `{ \tau^i H }`, where `i` ranges from 0 to `degree`.
    pub powers_of_h: Vec<G2Affine>,
    /// Pairing e(g, h) precomputed once during setup.
    pub e_gh: Gt,
}

#[derive(Debug)]
pub enum Error {
    /// The degree provided in setup was too small; degree 0 polynomials
    /// are not supported.
    DegreeIsZero,

    /// The degree of the polynomial passed to `commit` or `open`
    /// was too large.
    TooManyCoefficients {
        /// The number of coefficients in the polynomial.
        num_coefficients: usize,
        /// The maximum number of powers provided in `Powers`.
        num_powers: usize,
    },
}

impl KZG10 {
    pub fn setup(max_degree: usize, tau: Scalar) -> Result<PowersOfTau, Error> {
        if max_degree < 1 {
            return Err(Error::DegreeIsZero);
        }

        let g = G1Projective::generator();
        let h = G2Projective::generator();

        let mut powers_of_tau = vec![Scalar::ONE];

        let mut cur = tau;
        for _ in 0..max_degree {
            powers_of_tau.push(cur);
            cur *= tau;
        }
        
        // Compute powers of g: [g, g^tau, g^{tau^2}, ..., g^{tau^{max_degree}}]
        let powers_of_g_proj: Vec<G1Projective> =
            powers_of_tau.par_iter().map(|&power| g * power).collect();
        let mut powers_of_g = vec![G1Affine::default(); max_degree + 1];
        G1Projective::batch_normalize(&powers_of_g_proj, &mut powers_of_g);

        // Compute powers of h: [h, h^tau, h^{tau^2}, ..., h^{tau^{max_degree}}]
        let powers_of_h_proj: Vec<G2Projective> =
            powers_of_tau.par_iter().map(|&power| h * power).collect();
        let mut powers_of_h = vec![G2Affine::default(); max_degree + 1];
        G2Projective::batch_normalize(&powers_of_h_proj, &mut powers_of_h);

        // Precompute pairing e(g, h)
        let e_gh = blstrs::Bls12::pairing(&G1Affine::from(g), &G2Affine::from(h));

        let pp = PowersOfTau {
            powers_of_g,
            powers_of_h,
            e_gh,
        };

        Ok(pp)
    }

    pub fn commit_g1(
        params: &PowersOfTau,
        polynomial: &DensePolynomial,
    ) -> Result<G1Affine, Error> {
        let d = polynomial.degree();
        check_degree_is_too_large(d, params.powers_of_g.len())?;

        // MSM: sum of coeffs[i] * powers_of_g[i]
        let scalars = &polynomial.coeffs[..=d];
        let bases: Vec<G1Projective> = params.powers_of_g[..=d]
            .iter()
            .map(G1Projective::from)
            .collect();
        let commitment = G1Projective::multi_exp(&bases, scalars);

        Ok(commitment.to_affine())
    }

    pub fn commit_g2(
        params: &PowersOfTau,
        polynomial: &DensePolynomial,
    ) -> Result<G2Affine, Error> {
        let d = polynomial.degree();
        check_degree_is_too_large(d, params.powers_of_h.len())?;

        // MSM: sum of coeffs[i] * powers_of_h[i]
        let scalars = &polynomial.coeffs[..=d];
        let bases: Vec<G2Projective> = params.powers_of_h[..=d]
            .iter()
            .map(G2Projective::from)
            .collect();
        let commitment = G2Projective::multi_exp(&bases, scalars);

        Ok(commitment.to_affine())
    }

    pub fn compute_opening_proof(
        params: &PowersOfTau,
        polynomial: &DensePolynomial,
        point: &Scalar,
    ) -> Result<G1Affine, Error> {
        let eval = polynomial.evaluate(point);
        let eval_as_poly = DensePolynomial::from_coefficients_vec(vec![eval]);
        let numerator = polynomial.clone() - eval_as_poly;
        let divisor = DensePolynomial::from_coefficients_vec(vec![-*point, Scalar::ONE]);
        let witness_polynomial = &numerator / &divisor;

        Self::commit_g1(params, &witness_polynomial)
    }
}

fn check_degree_is_too_large(degree: usize, num_powers: usize) -> Result<(), Error> {
    let num_coefficients = degree + 1;
    if num_coefficients > num_powers {
        Err(Error::TooManyCoefficients {
            num_coefficients,
            num_powers,
        })
    } else {
        Ok(())
    }
}
