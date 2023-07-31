// use crate::utils::{lagrange_coefficients, transpose};
use ark_ec::{pairing::Pairing, Group};
use ark_poly::DenseUVPolynomial;
use ark_poly::{domain::EvaluationDomain, univariate::DensePolynomial, Radix2EvaluationDomain};
use ark_serialize::*;
use ark_std::{rand::RngCore, UniformRand, Zero};
use std::ops::{Mul, Sub};

use crate::encryption::Ciphertext;
use crate::kzg::{UniversalParams, KZG10};
use crate::utils::lagrange_poly;

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct SecretKey<E: Pairing> {
    sk: E::ScalarField,
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
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
    pub ask: E::G1,
    pub z_g2: E::G2,
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

    pub fn get_pk(&self, id: usize, params: &UniversalParams<E>, n: usize) -> PublicKey<E> {
        let domain = Radix2EvaluationDomain::<E::ScalarField>::new(n).unwrap();

        let li = lagrange_poly(n, id);

        let mut sk_li_by_z = vec![];
        for j in 0..n {
            let num: DensePolynomial<E::ScalarField>;
            if id == j {
                num = li.clone().mul(&li).sub(&li);
            } else {
                //cross-terms
                let l_j = lagrange_poly(n, j);
                num = l_j.mul(&li);
            }

            let f = num.divide_by_vanishing_poly(domain).unwrap().0;
            let sk_times_f = &f * self.sk;

            let com = KZG10::commit_g1(&params, &sk_times_f)
                .expect("commitment failed")
                .into();

            sk_li_by_z.push(com);
        }

        let f = DensePolynomial::from_coefficients_vec(li.coeffs[1..].to_vec());
        let sk_times_f = &f * self.sk;
        let sk_li_by_tau = KZG10::commit_g1(&params, &sk_times_f)
            .expect("commitment failed")
            .into();

        let mut f = &li * self.sk;
        let sk_li = KZG10::commit_g1(&params, &f)
            .expect("commitment failed")
            .into();

        f.coeffs[0] = E::ScalarField::zero();
        let sk_li_minus0 = KZG10::commit_g1(&params, &f)
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
        let domain = Radix2EvaluationDomain::<E::ScalarField>::new(pk.len()).unwrap();

        // todo: replace this with efficient z_g2.
        let z = domain.vanishing_polynomial().clone();
        let z_g2 = KZG10::<E, DensePolynomial<E::ScalarField>>::commit_g2(&params, &z.into())
            .expect("commitment failed")
            .into();

        // gather sk_li from all public keys
        let mut ask = E::G1::zero();
        for i in 0..pk.len() {
            ask += pk[i].sk_li;
        }

        AggregateKey { pk, ask, z_g2 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type E = ark_bls12_381::Bls12_381;
    type UniPoly381 = DensePolynomial<<E as Pairing>::ScalarField>;
    // type F = <E as Pairing>::ScalarField;
    // type G1 = <E as Pairing>::G1;
    // type G2 = <E as Pairing>::G2;

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

// fn party_i_setup_material(
//     params: &UniversalParams<Bls12_381>,
//     n: usize,
//     i: usize,
//     sk_i: &F) -> (G1, G2, Vec<G1>, G1) {
//     //let us compute the q1 term
//     let l_i_of_x = utils::lagrange_poly(n, i);
//     let z_of_x = utils::compute_vanishing_poly(n);

//     let mut q1_material = vec![];
//     //let us compute the cross terms of q1
//     for j in 0..n {
//         let num: DensePolynomial<F>;// = compute_constant_poly(&F::from(0));
//         if i == j {
//             num = l_i_of_x.clone().mul(&l_i_of_x).sub(&l_i_of_x);
//         } else { //cross-terms
//             let l_j_of_x = utils::lagrange_poly(n, j);
//             num = l_j_of_x.mul(&l_i_of_x);
//         }
//         let f = num.div(&z_of_x);
//         let sk_times_f = utils::poly_eval_mult_c(&f, &sk_i);

//         let com = KZG::commit_g1(&params, &sk_times_f)
//             .expect("commitment failed");

//         q1_material.push(com);
//     }

//     let x_monomial = utils::compute_x_monomial();
//     let l_i_of_0 = l_i_of_x.evaluate(&F::from(0));
//     let l_i_of_0_poly = utils::compute_constant_poly(&l_i_of_0);
//     let num = l_i_of_x.sub(&l_i_of_0_poly);
//     let den = x_monomial.clone();
//     let f = num.div(&den);
//     let sk_times_f = utils::poly_eval_mult_c(&f, &sk_i);
//     let q2_com = KZG::commit_g1(&params, &sk_times_f).expect("commitment failed");

//     //release my public key
//     let sk_as_poly = utils::compute_constant_poly(&sk_i);
//     let pk = KZG::commit_g1(&params, &sk_as_poly).expect("commitment failed");

//     let sk_times_l_i_of_x = utils::poly_eval_mult_c(&l_i_of_x, &sk_i);
//     let com_sk_l_i = KZG::commit_g2(&params, &sk_times_l_i_of_x).expect("commitment failed");

//     (pk, com_sk_l_i, q1_material, q2_com)
// }
