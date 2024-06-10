//adapted from https://github.com/arkworks-rs/poly-commit/blob/master/src/kzg10/mod.rs
#![allow(dead_code)]
#![allow(unused_imports)]

use std::cmp::Ordering;
use std::collections::HashMap;

use ark_ec::{pairing::Pairing, CurveGroup, Group};
use ark_ec::{scalar_mul::fixed_base::FixedBase, VariableBaseMSM};
use ark_ff::{One, PrimeField, UniformRand, Zero};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{
    DenseUVPolynomial, EvaluationDomain as _, Evaluations, Polynomial as _, Radix2EvaluationDomain,
};
use ark_std::{format, marker::PhantomData, ops::*, vec};

use ark_std::rand::RngCore;
use rayon::iter::IntoParallelIterator as _;
use rayon::iter::ParallelIterator as _;

pub struct KZG10<E: Pairing> {
    _engine: PhantomData<E>,
}

pub struct UniversalParams<E: Pairing> {
    /// Group elements of the form `{ \beta^i G }`, where `i` ranges from 0 to `degree`.
    pub powers_of_g: Vec<E::G1Affine>,
    /// Group elements of the form `{ \beta^i H }`, where `i` ranges from 0 to `degree`.
    pub powers_of_h: Vec<E::G2Affine>,

    // Lagrange polynomials
    pub l_i: Vec<DensePolynomial<E::ScalarField>>,
    // Precomputed commitments
    pub li_by_z: HashMap<(usize, usize), E::G1>,
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

impl<E> KZG10<E>
where
    E: Pairing,
{
    pub fn setup<R: RngCore>(max_degree: usize, rng: &mut R) -> Result<UniversalParams<E>, Error> {
        if max_degree < 1 {
            return Err(Error::DegreeIsZero);
        }

        //let setup_time = start_timer!(|| format!("KZG10::Setup with degree {}", max_degree));
        let beta = E::ScalarField::rand(rng);
        let g = E::G1::generator();
        let h = E::G2::generator();

        let mut powers_of_beta = vec![E::ScalarField::one()];

        let mut cur = beta;
        for _ in 0..max_degree {
            powers_of_beta.push(cur);
            cur *= &beta;
        }

        let window_size = FixedBase::get_mul_window_size(max_degree + 1);
        let scalar_bits = E::ScalarField::MODULUS_BIT_SIZE as usize;

        let g_table = FixedBase::get_window_table(scalar_bits, window_size, g);
        let powers_of_g =
            FixedBase::msm::<E::G1>(scalar_bits, window_size, &g_table, &powers_of_beta);

        let h_table = FixedBase::get_window_table(scalar_bits, window_size, h);
        let powers_of_h =
            FixedBase::msm::<E::G2>(scalar_bits, window_size, &h_table, &powers_of_beta);

        let powers_of_g = E::G1::normalize_batch(&powers_of_g);
        let powers_of_h = E::G2::normalize_batch(&powers_of_h);

        let l_i = (0..max_degree)
            .map(|i| Self::lagrange_poly(max_degree, i))
            .collect::<Vec<_>>();
        let mut pp = UniversalParams {
            powers_of_g,
            powers_of_h,
            li_by_z: HashMap::new(),
            l_i,
        };
        let domain = Radix2EvaluationDomain::<E::ScalarField>::new(max_degree).unwrap();

        pp.li_by_z = (0..max_degree)
            .into_par_iter()
            .map(|i| (i..max_degree).into_par_iter().map(move |j| (i, j)))
            .flatten()
            .filter_map(|(i, j)| match i.cmp(&j) {
                Ordering::Less => Some((i, j, pp.l_i[j].mul(&pp.l_i[i]))),
                Ordering::Equal => Some((i, j, pp.l_i[i].mul(&pp.l_i[i]).sub(&pp.l_i[i]))),
                Ordering::Greater => None,
            })
            .map(|(i, j, numerator)| (i, j, numerator.divide_by_vanishing_poly(domain).unwrap().0))
            .map(|(i, j, f)| ((i, j), Self::commit_g1(&pp, &f).unwrap().into()))
            .collect::<HashMap<(usize, usize), E::G1>>();
        //end_timer!(setup_time);
        Ok(pp)
    }

    pub fn commit_g1(
        params: &UniversalParams<E>,
        polynomial: &DensePolynomial<E::ScalarField>,
    ) -> Result<E::G1Affine, Error> {
        let d = polynomial.degree();
        check_degree_is_too_large(d, params.powers_of_g.len())?;

        let plain_coeffs = convert_to_bigints(polynomial.coeffs());

        let powers_of_g = &params.powers_of_g[..=d].to_vec();
        //let msm_time = start_timer!(|| "MSM to compute commitment to plaintext poly");
        let commitment = <E::G1 as VariableBaseMSM>::msm_bigint(&powers_of_g[..], &plain_coeffs);
        //end_timer!(msm_time);
        Ok(commitment.into_affine())
    }

    pub fn commit_g2(
        params: &UniversalParams<E>,
        polynomial: &DensePolynomial<E::ScalarField>,
    ) -> Result<E::G2Affine, Error> {
        let d = polynomial.degree();
        check_degree_is_too_large(d, params.powers_of_h.len())?;

        let plain_coeffs = convert_to_bigints(polynomial.coeffs());

        let powers_of_h = &params.powers_of_h[..=d].to_vec();
        //let msm_time = start_timer!(|| "MSM to compute commitment to plaintext poly");
        let commitment = <E::G2 as VariableBaseMSM>::msm_bigint(&powers_of_h[..], &plain_coeffs);
        //end_timer!(msm_time);

        Ok(commitment.into_affine())
    }

    // 1 at omega^i and 0 elsewhere on domain {omega^i}_{i \in [n]}
    pub fn lagrange_poly(n: usize, i: usize) -> DensePolynomial<E::ScalarField> {
        debug_assert!(i < n);
        //todo: check n is a power of 2
        let mut evals = vec![];
        for j in 0..n {
            let l_of_x: u64 = if i == j { 1 } else { 0 };
            evals.push(E::ScalarField::from(l_of_x));
        }

        //powers of nth root of unity
        let domain = Radix2EvaluationDomain::<E::ScalarField>::new(n).unwrap();
        let eval_form = Evaluations::from_vec_and_domain(evals, domain);
        //interpolated polynomial over the n points
        eval_form.interpolate()
    }

    pub fn compute_opening_proof(
        params: &UniversalParams<E>,
        polynomial: &DensePolynomial<E::ScalarField>,
        point: &E::ScalarField,
    ) -> Result<E::G1Affine, Error> {
        let eval = polynomial.evaluate(point);

        let eval_as_poly =
            <DensePolynomial<E::ScalarField> as DenseUVPolynomial<E::ScalarField>>::from_coefficients_vec(vec![
                eval,
            ]);
        let numerator = polynomial.clone().sub(&eval_as_poly);
        let divisor = DenseUVPolynomial::<E::ScalarField>::from_coefficients_vec(vec![
            E::ScalarField::zero() - point,
            E::ScalarField::one(),
        ]);
        let witness_polynomial = numerator.div(&divisor);
        Self::commit_g1(params, &witness_polynomial)
    }
}

fn skip_leading_zeros_and_convert_to_bigints<F: PrimeField, P: DenseUVPolynomial<F>>(
    p: &P,
) -> (usize, Vec<F::BigInt>) {
    let mut num_leading_zeros = 0;
    while num_leading_zeros < p.coeffs().len() && p.coeffs()[num_leading_zeros].is_zero() {
        num_leading_zeros += 1;
    }
    let coeffs = convert_to_bigints(&p.coeffs()[num_leading_zeros..]);
    (num_leading_zeros, coeffs)
}

pub fn convert_to_bigints<F: PrimeField>(p: &[F]) -> Vec<F::BigInt> {
    ark_std::cfg_iter!(p)
        .map(|s| s.into_bigint())
        .collect::<Vec<_>>()
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
