use ark_ec::pairing::PairingOutput;
// use crate::utils::{lagrange_coefficients, transpose};
use crate::encryption::Ciphertext;
use crate::kzg::{PowersOfTau, KZG10};
use crate::utils::lagrange_poly;
use ark_ec::{pairing::Pairing, PrimeGroup};
use ark_ff::Field;
use ark_poly::{
    domain::EvaluationDomain, univariate::DensePolynomial, DenseUVPolynomial, Polynomial,
    Radix2EvaluationDomain,
};
use ark_serialize::*;
use ark_std::{rand::RngCore, One, UniformRand, Zero};
use rayon::prelude::*;
use std::ops::{Mul, Sub};

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct LagrangePowers<E: Pairing> {
    pub li: Vec<E::G1>,
    pub li_minus0: Vec<E::G1>,
    pub li_x: Vec<E::G1>,
    pub li_lj_z: Vec<Vec<E::G1>>,
}

impl<E: Pairing> LagrangePowers<E> {
    pub fn new(tau: E::ScalarField, n: usize) -> Self {
        let mut li_evals: Vec<E::ScalarField> = vec![E::ScalarField::zero(); n];
        let mut li_evals_minus0: Vec<E::ScalarField> = vec![E::ScalarField::zero(); n];
        let mut li_evals_x: Vec<E::ScalarField> = vec![E::ScalarField::zero(); n];
        let tau_inv = tau.inverse().unwrap();
        for i in 0..n {
            let li = lagrange_poly(n, i);
            li_evals[i] = li.evaluate(&tau);

            li_evals_minus0[i] = li_evals[i] - li.coeffs[0];

            li_evals_x[i] = li_evals_minus0[i] * tau_inv;
        }

        let z_eval = tau.pow(&[n as u64]) - E::ScalarField::one();
        let z_eval_inv = z_eval.inverse().unwrap();

        let mut li = vec![E::G1::zero(); n];
        for i in 0..n {
            li[i] = E::G1::generator() * li_evals[i];
        }

        let mut li_minus0 = vec![E::G1::zero(); n];
        li_minus0.par_iter_mut().enumerate().for_each(|(i, elem)| {
            *elem = E::G1::generator() * li_evals_minus0[i];
        });

        let mut li_x = vec![E::G1::zero(); n];
        li_x.par_iter_mut().enumerate().for_each(|(i, elem)| {
            *elem = E::G1::generator() * li_evals_x[i];
        });

        let mut li_lj_z = vec![vec![E::G1::zero(); n]; n];
        li_lj_z.par_iter_mut().enumerate().for_each(|(i, row)| {
            row.par_iter_mut().enumerate().for_each(|(j, elem)| {
                *elem = if i == j {
                    E::G1::generator() * ((li_evals[i] * li_evals[i] - li_evals[i]) * z_eval_inv)
                } else {
                    E::G1::generator() * (li_evals[i] * li_evals[j] * z_eval_inv)
                }
            });
        });

        LagrangePowers {
            li,
            li_minus0,
            li_x,
            li_lj_z,
        }
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct SecretKey<E: Pairing> {
    sk: E::ScalarField,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Default, Debug)]
pub struct PublicKey<E: Pairing> {
    pub id: usize,
    pub bls_pk: E::G1,          //BLS pk
    pub sk_li: E::G1,           //hint
    pub sk_li_minus0: E::G1,    //hint
    pub sk_li_lj_z: Vec<E::G1>, //hint
    pub sk_li_x: E::G1,         //hint
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct AggregateKey<E: Pairing> {
    pub pk: Vec<PublicKey<E>>,
    pub agg_sk_li_lj_z: Vec<E::G1>,
    pub ask: E::G1,
    pub z_g2: E::G2,

    //preprocessed values
    pub h_minus1: E::G2,
    pub e_gh: PairingOutput<E>,
}

impl<E: Pairing> PublicKey<E> {
    pub fn new(
        id: usize,
        bls_pk: E::G1,
        sk_li: E::G1,
        sk_li_minus0: E::G1,
        sk_li_lj_z: Vec<E::G1>,
        sk_li_x: E::G1,
    ) -> Self {
        PublicKey {
            id,
            bls_pk,
            sk_li,
            sk_li_minus0,
            sk_li_lj_z,
            sk_li_x,
        }
    }
}

impl<E: Pairing> SecretKey<E> {
    pub fn new<R: RngCore>(rng: &mut R) -> Self {
        SecretKey {
            sk: E::ScalarField::rand(rng),
        }
    }

    pub fn nullify(&mut self) {
        self.sk = E::ScalarField::one()
    }

    pub fn get_pk(&self, id: usize, params: &PowersOfTau<E>, n: usize) -> PublicKey<E> {
        // TODO: This runs in quadratic time because we are not preprocessing the Li's
        // Fix this.
        let domain = Radix2EvaluationDomain::<E::ScalarField>::new(n).unwrap();

        let li = lagrange_poly(n, id);

        let mut sk_li_lj_z = vec![];
        for j in 0..n {
            let num = if id == j {
                li.clone().mul(&li).sub(&li)
            } else {
                //cross-terms
                let l_j = lagrange_poly(n, j);
                l_j.mul(&li)
            };

            let f = num.divide_by_vanishing_poly(domain).0;
            let sk_times_f = &f * self.sk;

            let com = KZG10::commit_g1(params, &sk_times_f)
                .expect("commitment failed")
                .into();

            sk_li_lj_z.push(com);
        }

        let f = DensePolynomial::from_coefficients_vec(li.coeffs[1..].to_vec());
        let sk_times_f = &f * self.sk;
        let sk_li_x = KZG10::commit_g1(params, &sk_times_f)
            .expect("commitment failed")
            .into();

        let mut f = &li * self.sk;
        let sk_li = KZG10::commit_g1(params, &f)
            .expect("commitment failed")
            .into();

        f.coeffs[0] = E::ScalarField::zero();
        let sk_li_minus0 = KZG10::commit_g1(params, &f)
            .expect("commitment failed")
            .into();

        PublicKey {
            id,
            bls_pk: E::G1::generator() * self.sk,
            sk_li,
            sk_li_minus0,
            sk_li_lj_z,
            sk_li_x,
        }
    }

    pub fn lagrange_get_pk(&self, id: usize, params: &LagrangePowers<E>, n: usize) -> PublicKey<E> {
        let mut sk_li_lj_z = vec![];

        let sk_li = params.li[id] * self.sk;

        let sk_li_minus0 = params.li_minus0[id] * self.sk;

        let sk_li_x = params.li_x[id] * self.sk;

        for j in 0..n {
            sk_li_lj_z.push(params.li_lj_z[id][j] * self.sk);
        }

        PublicKey {
            id,
            bls_pk: E::G1::generator() * self.sk,
            sk_li,
            sk_li_minus0,
            sk_li_lj_z,
            sk_li_x,
        }
    }

    pub fn partial_decryption(&self, ct: &Ciphertext<E>) -> E::G2 {
        ct.gamma_g2 * self.sk // kind of a bls signature on gamma_g2
    }
}

impl<E: Pairing> AggregateKey<E> {
    pub fn new(pk: Vec<PublicKey<E>>, params: &PowersOfTau<E>) -> Self {
        let n = pk.len();
        let h_minus1 = params.powers_of_h[0] * (-E::ScalarField::one());
        let z_g2 = params.powers_of_h[n] + h_minus1;

        // gather sk_li from all public keys
        let mut ask = E::G1::zero();
        for pki in pk.iter() {
            ask += pki.sk_li;
        }

        let mut agg_sk_li_lj_z = vec![];
        for i in 0..n {
            let mut agg_sk_li_lj_zi = E::G1::zero();
            for pkj in pk.iter() {
                agg_sk_li_lj_zi += pkj.sk_li_lj_z[i];
            }
            agg_sk_li_lj_z.push(agg_sk_li_lj_zi);
        }

        AggregateKey {
            pk,
            agg_sk_li_lj_z,
            ask,
            z_g2,
            h_minus1,
            e_gh: E::pairing(params.powers_of_g[0], params.powers_of_h[0]),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type E = ark_bls12_381::Bls12_381;
    type Fr = <E as Pairing>::ScalarField;
    type UniPoly381 = DensePolynomial<<E as Pairing>::ScalarField>;

    #[test]
    fn test_setup() {
        let mut rng = ark_std::test_rng();
        let n = 16;
        let tau = Fr::rand(&mut rng);
        let params = KZG10::<E, UniPoly381>::setup(n, tau.clone()).unwrap();
        let lagrange_params = LagrangePowers::<E>::new(tau, n);

        let mut sk: Vec<SecretKey<E>> = Vec::new();
        let mut pk: Vec<PublicKey<E>> = Vec::new();
        let mut lagrange_pk: Vec<PublicKey<E>> = Vec::new();

        for i in 0..n {
            sk.push(SecretKey::<E>::new(&mut rng));
            pk.push(sk[i].get_pk(i, &params, n));
            lagrange_pk.push(sk[i].lagrange_get_pk(i, &lagrange_params, n));

            assert_eq!(pk[i].sk_li, lagrange_pk[i].sk_li);
            assert_eq!(pk[i].sk_li_minus0, lagrange_pk[i].sk_li_minus0);
            assert_eq!(pk[i].sk_li_x, lagrange_pk[i].sk_li_x); //computed incorrectly go fix it
            assert_eq!(pk[i].sk_li_lj_z, lagrange_pk[i].sk_li_lj_z);
        }

        let _ak = AggregateKey::<E>::new(pk, &params);
    }
}
