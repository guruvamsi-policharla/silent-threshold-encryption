use ark_ec::pairing::PairingOutput;
// use crate::utils::{lagrange_coefficients, transpose};
use ark_ec::{pairing::Pairing, Group};
use ark_poly::DenseUVPolynomial;
use ark_poly::{domain::EvaluationDomain, univariate::DensePolynomial, Radix2EvaluationDomain};
use ark_serialize::*;
use ark_std::{rand::RngCore, One, UniformRand, Zero};
use std::ops::{Mul, Sub};

use crate::encryption::Ciphertext;
use crate::kzg::{UniversalParams, KZG10};
use crate::utils::lagrange_poly;

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct SecretKey<E: Pairing> {
    sk: E::ScalarField,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct PublicKey<E: Pairing> {
    pub id: usize,
    pub bls_pk: E::G1,          //BLS pk
    pub sk_li: E::G1,           //hint
    pub sk_li_minus0: E::G1,    //hint
    pub sk_li_by_z: Vec<E::G1>, //hint
    pub sk_li_by_tau: E::G1,    //hint
}

pub struct AggregateKey<E: Pairing> {
    pub pk: Vec<PublicKey<E>>,
    pub agg_sk_li_by_z: Vec<E::G1>,
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
        sk_li_by_z: Vec<E::G1>,
        sk_li_by_tau: E::G1,
    ) -> Self {
        PublicKey {
            id,
            bls_pk,
            sk_li,
            sk_li_minus0,
            sk_li_by_z,
            sk_li_by_tau,
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

    pub fn get_pk(&self, id: usize, params: &UniversalParams<E>, n: usize) -> PublicKey<E> {
        let domain = Radix2EvaluationDomain::<E::ScalarField>::new(n).unwrap();

        let li = lagrange_poly(n, id);

        let mut sk_li_by_z = vec![];
        for j in 0..n {
            let num = if id == j {
                li.clone().mul(&li).sub(&li)
            } else {
                //cross-terms
                let l_j = lagrange_poly(n, j);
                l_j.mul(&li)
            };

            let f = num.divide_by_vanishing_poly(domain).unwrap().0;
            let sk_times_f = &f * self.sk;

            let com = KZG10::commit_g1(params, &sk_times_f)
                .expect("commitment failed")
                .into();

            sk_li_by_z.push(com);
        }

        let f = DensePolynomial::from_coefficients_vec(li.coeffs[1..].to_vec());
        let sk_times_f = &f * self.sk;
        let sk_li_by_tau = KZG10::commit_g1(params, &sk_times_f)
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
            sk_li_by_z,
            sk_li_by_tau,
        }
    }

    pub fn partial_decryption(&self, ct: &Ciphertext<E>) -> E::G2 {
        ct.gamma_g2 * self.sk // kind of a bls signature on gamma_g2
    }
}

impl<E: Pairing> AggregateKey<E> {
    pub fn new(pk: Vec<PublicKey<E>>, params: &UniversalParams<E>) -> Self {
        let n = pk.len();
        let h_minus1 = params.powers_of_h[0] * (-E::ScalarField::one());
        let z_g2 = params.powers_of_h[n] + h_minus1;

        // gather sk_li from all public keys
        let mut ask = E::G1::zero();
        for pki in pk.iter() {
            ask += pki.sk_li;
        }

        let mut agg_sk_li_by_z = vec![];
        for i in 0..n {
            let mut agg_sk_li_by_zi = E::G1::zero();
            for pkj in pk.iter() {
                agg_sk_li_by_zi += pkj.sk_li_by_z[i];
            }
            agg_sk_li_by_z.push(agg_sk_li_by_zi);
        }

        AggregateKey {
            pk,
            agg_sk_li_by_z,
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
    type UniPoly381 = DensePolynomial<<E as Pairing>::ScalarField>;

    #[test]
    fn test_setup() {
        let mut rng = ark_std::test_rng();
        let n = 4;
        let params = KZG10::<E, UniPoly381>::setup(n, &mut rng).unwrap();

        let mut sk: Vec<SecretKey<E>> = Vec::new();
        let mut pk: Vec<PublicKey<E>> = Vec::new();

        for i in 0..n {
            sk.push(SecretKey::<E>::new(&mut rng));
            pk.push(sk[i].get_pk(0, &params, n))
        }

        let _ak = AggregateKey::<E>::new(pk, &params);
    }
}
