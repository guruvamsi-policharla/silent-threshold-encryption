//adapted from https://github.com/arkworks-rs/poly-commit/blob/master/src/kzg10/mod.rs
#![allow(dead_code)]
#![allow(unused_imports)]

use blstrs::{G1Affine, G1Projective, G2Affine, G2Projective, Scalar};
use ff::{Field, PrimeField};
use group::{Curve, Group};
use serde::{Deserialize, Serialize};
use std::ops::*;

use crate::polynomial::DensePolynomial;

pub struct KZG10;

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct PowersOfTau {
    /// Group elements of the form `{ \tau^i G }`, where `i` ranges from 0 to `degree`.
    pub powers_of_g: Vec<G1Affine>,
    /// Group elements of the form `{ \tau^i H }`, where `i` ranges from 0 to `degree`.
    pub powers_of_h: Vec<G2Affine>,
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
        for _ in 0..=max_degree {
            powers_of_tau.push(cur);
            cur *= tau;
        }

        // Compute powers of g: [g, g^tau, g^{tau^2}, ..., g^{tau^{max_degree}}]
        let powers_of_g: Vec<G1Affine> = powers_of_tau[0..=max_degree]
            .iter()
            .map(|&power| (g * power).to_affine())
            .collect();

        // Compute powers of h: [h, h^tau, h^{tau^2}, ..., h^{tau^{max_degree}}]
        let powers_of_h: Vec<G2Affine> = powers_of_tau[0..=max_degree]
            .iter()
            .map(|&power| (h * power).to_affine())
            .collect();

        let pp = PowersOfTau {
            powers_of_g,
            powers_of_h,
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
        let commitment = polynomial
            .coeffs
            .iter()
            .zip(&params.powers_of_g[..=d])
            .fold(G1Projective::identity(), |acc, (coeff, base)| {
                acc + base * coeff
            });

        Ok(commitment.to_affine())
    }

    pub fn commit_g2(
        params: &PowersOfTau,
        polynomial: &DensePolynomial,
    ) -> Result<G2Affine, Error> {
        let d = polynomial.degree();
        check_degree_is_too_large(d, params.powers_of_h.len())?;

        // MSM: sum of coeffs[i] * powers_of_h[i]
        let commitment = polynomial
            .coeffs
            .iter()
            .zip(&params.powers_of_h[..=d])
            .fold(G2Projective::identity(), |acc, (coeff, base)| {
                acc + base * coeff
            });

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
