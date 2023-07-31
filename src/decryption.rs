use std::ops::Div;

use ark_ec::pairing::Pairing;
use ark_ff::Field;
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain,
};
use ark_std::{One, Zero};

use crate::{
    encryption::Ciphertext,
    kzg::{UniversalParams, KZG10},
    setup::{AggregateKey, SecretKey},
    utils::interp_mostly_zero,
};

pub fn agg_dec<E: Pairing>(
    partial_decryptions: &Vec<E::G2>,
    ct: &Ciphertext<E>,
    selector: &Vec<bool>,
    agg_key: &AggregateKey<E>,
    params: &UniversalParams<E>,
) {
    // using iter create a DensePolynomial b from selector where b[i] = E::ScalarField::one() if selector[i] = true
    // and b[i] = E::ScalarField::zero() if selector[i] = false
    let domain = Radix2EvaluationDomain::<E::ScalarField>::new(selector.len()).unwrap();
    let domain_elements: Vec<E::ScalarField> = domain.elements().collect();

    // select the points that have selector[i] = true
    let mut points: Vec<E::ScalarField> = Vec::new();
    let mut parties: Vec<usize> = Vec::new(); // parties indexed from 0..n-1
    selector.iter().enumerate().for_each(|(i, x)| {
        if *x {
            points.push(domain_elements[i]);
            parties.push(i);
        }
    });

    let b = interp_mostly_zero(E::ScalarField::one(), &points);
    let b_evals = domain.fft(&b.coeffs);

    // commit to b in g2
    let b_g2 = KZG10::<E, DensePolynomial<E::ScalarField>>::commit_g2(&params, &b)
        .unwrap()
        .into();

    // accumulate the public key
    // count number of true in selector
    let mut v: usize = 0;
    selector.iter().for_each(|x| {
        if *x {
            v += 1;
        }
    });
    let mut v = E::ScalarField::from(v as u32);
    v.inverse_in_place();

    // compute the aggregate public key
    let mut apk = E::G1::zero();
    for i in 0..agg_key.pk.len() {
        if selector[i] {
            apk += agg_key.pk[i].bls_pk;
        }
    }
    apk *= v;

    // compute Qx, Qhatx and Qz
    let mut qx: E::G1 = E::G1::zero();
    for &i in &parties {
        qx += agg_key.pk[i].sk_li_by_tau * b_evals[i];
    }

    let mut qz: E::G1 = E::G1::zero();
    for &i in &parties {
        let mut qz_i = E::G1::zero();
        for &j in &parties {
            qz_i += agg_key.pk[i].sk_li_by_z[j];
        }
        qz += qz_i * b_evals[i];
    }

    let mut qhatx: E::G1 = E::G1::zero();
    for &i in &parties {
        qhatx += agg_key.pk[i].sk_li_minus0 * b_evals[i];
    }

    // compute sigma = (\sum partial_decryptions[i])/v for i in parties
    let mut sigma: E::G1 = E::G1::zero();
    for &i in &parties {
        sigma += partial_decryptions[i];
    }
    sigma *= v;

    // bhat = x^t * b
    // insert t 0s at the beginning of bhat.coeffs
    let mut bhat_coeffs = vec![E::ScalarField::zero(); ct.t];
    bhat_coeffs.append(&mut b.coeffs.clone());
    let bhat = KZG10::<E, DensePolynomial<E::ScalarField>>::commit_g1(
        &params,
        &DensePolynomial::from_coefficients_vec(bhat_coeffs),
    )
    .unwrap()
    .into();

    // q0 = (b-1)/(x-domain_elements[0])
    let mut bminus1 = b.clone();
    bminus1.coeffs[0] -= E::ScalarField::one();

    let q0 = b.div(&DensePolynomial::from_coefficients_vec(vec![
        -domain_elements[0],
        E::ScalarField::one(),
    ]));

    let q0_g1 = KZG10::<E, DensePolynomial<E::ScalarField>>::commit_g1(&params, &q0)
        .unwrap()
        .into();

    let w1 = [apk, qz, qx, qhatx, bhat, q0_g1];
    let w2 = [sigma, b_g2];
}
