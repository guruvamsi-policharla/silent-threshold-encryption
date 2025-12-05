use blstrs::{G1Affine, G1Projective, G2Affine, G2Projective, Gt, Scalar};
use ff::Field;
use group::Group;
use serde::{Deserialize, Serialize};

use crate::encryption::Ciphertext;
use crate::kzg::{PowersOfTau, KZG10};
use crate::polynomial::{DensePolynomial, Radix2EvaluationDomain};
use crate::utils::lagrange_poly;
use rayon::prelude::*;

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct LagrangePowers {
    pub li: Vec<G1Projective>,
    pub li_minus0: Vec<G1Projective>,
    pub li_x: Vec<G1Projective>,
    pub li_lj_z: Vec<Vec<G1Projective>>,
}

impl LagrangePowers {
    pub fn new(tau: Scalar, n: usize) -> Self {
        let mut li_evals: Vec<Scalar> = vec![Scalar::ZERO; n];
        let mut li_evals_minus0: Vec<Scalar> = vec![Scalar::ZERO; n];
        let mut li_evals_x: Vec<Scalar> = vec![Scalar::ZERO; n];
        let tau_inv = tau.invert().unwrap();
        for i in 0..n {
            let li = lagrange_poly(n, i);
            li_evals[i] = li.evaluate(&tau);

            li_evals_minus0[i] = li_evals[i] - li.coeffs[0];

            li_evals_x[i] = li_evals_minus0[i] * tau_inv;
        }

        let z_eval = tau.pow(&[n as u64, 0, 0, 0]) - Scalar::ONE;
        let z_eval_inv = z_eval.invert().unwrap();

        let mut li = vec![G1Projective::identity(); n];
        for i in 0..n {
            li[i] = G1Projective::generator() * li_evals[i];
        }

        let mut li_minus0 = vec![G1Projective::identity(); n];
        li_minus0.par_iter_mut().enumerate().for_each(|(i, elem)| {
            *elem = G1Projective::generator() * li_evals_minus0[i];
        });

        let mut li_x = vec![G1Projective::identity(); n];
        li_x.par_iter_mut().enumerate().for_each(|(i, elem)| {
            *elem = G1Projective::generator() * li_evals_x[i];
        });

        let mut li_lj_z = vec![vec![G1Projective::identity(); n]; n];
        li_lj_z.par_iter_mut().enumerate().for_each(|(i, row)| {
            row.par_iter_mut().enumerate().for_each(|(j, elem)| {
                *elem = if i == j {
                    G1Projective::generator()
                        * ((li_evals[i] * li_evals[i] - li_evals[i]) * z_eval_inv)
                } else {
                    G1Projective::generator() * (li_evals[i] * li_evals[j] * z_eval_inv)
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

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct SecretKey {
    sk: Scalar,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct PublicKey {
    pub id: usize,
    pub bls_pk: G1Projective,          //BLS pk
    pub sk_li: G1Projective,           //hint
    pub sk_li_minus0: G1Projective,    //hint
    pub sk_li_lj_z: Vec<G1Projective>, //hint
    pub sk_li_x: G1Projective,         //hint
}

impl Default for PublicKey {
    fn default() -> Self {
        PublicKey {
            id: 0,
            bls_pk: G1Projective::identity(),
            sk_li: G1Projective::identity(),
            sk_li_minus0: G1Projective::identity(),
            sk_li_lj_z: vec![],
            sk_li_x: G1Projective::identity(),
        }
    }
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct AggregateKey {
    pub pk: Vec<PublicKey>,
    pub agg_sk_li_lj_z: Vec<G1Projective>,
    pub ask: G1Projective,
    pub z_g2: G2Projective,

    //preprocessed values
    pub h_minus1: G2Projective,
    pub e_gh: Gt,
}

impl PublicKey {
    pub fn new(
        id: usize,
        bls_pk: G1Projective,
        sk_li: G1Projective,
        sk_li_minus0: G1Projective,
        sk_li_lj_z: Vec<G1Projective>,
        sk_li_x: G1Projective,
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

impl SecretKey {
    pub fn new<R: rand_core::RngCore>(rng: &mut R) -> Self {
        SecretKey {
            sk: Scalar::random(rng),
        }
    }

    pub fn nullify(&mut self) {
        self.sk = Scalar::ONE
    }

    pub fn get_pk(&self, id: usize, params: &PowersOfTau, n: usize) -> PublicKey {
        // TODO: This runs in quadratic time because we are not preprocessing the Li's
        // Fix this.
        let domain = Radix2EvaluationDomain::new(n).unwrap();

        let li = lagrange_poly(n, id);

        let mut sk_li_lj_z = vec![];
        for j in 0..n {
            let num = if id == j {
                li.clone().naive_mul(&li) - li.clone()
            } else {
                //cross-terms
                let l_j = lagrange_poly(n, j);
                l_j.naive_mul(&li)
            };

            let f = num.divide_by_vanishing_poly(domain.clone()).0;
            let sk_times_f = &f * self.sk;

            let com = KZG10::commit_g1(params, &sk_times_f).expect("commitment failed");

            sk_li_lj_z.push(com.into());
        }

        let f = DensePolynomial::from_coefficients_vec(li.coeffs[1..].to_vec());
        let sk_times_f = &f * self.sk;
        let sk_li_x: G1Projective = KZG10::commit_g1(params, &sk_times_f)
            .expect("commitment failed")
            .into();

        let mut f = &li * self.sk;
        let sk_li: G1Projective = KZG10::commit_g1(params, &f)
            .expect("commitment failed")
            .into();

        f.coeffs[0] = Scalar::ZERO;
        let sk_li_minus0: G1Projective = KZG10::commit_g1(params, &f)
            .expect("commitment failed")
            .into();

        PublicKey {
            id,
            bls_pk: G1Projective::generator() * self.sk,
            sk_li,
            sk_li_minus0,
            sk_li_lj_z,
            sk_li_x,
        }
    }

    pub fn lagrange_get_pk(&self, id: usize, params: &LagrangePowers, n: usize) -> PublicKey {
        let mut sk_li_lj_z = vec![];

        let sk_li = params.li[id] * self.sk;

        let sk_li_minus0 = params.li_minus0[id] * self.sk;

        let sk_li_x = params.li_x[id] * self.sk;

        for j in 0..n {
            sk_li_lj_z.push(params.li_lj_z[id][j] * self.sk);
        }

        PublicKey {
            id,
            bls_pk: G1Projective::generator() * self.sk,
            sk_li,
            sk_li_minus0,
            sk_li_lj_z,
            sk_li_x,
        }
    }

    pub fn partial_decryption(&self, ct: &Ciphertext) -> G2Projective {
        ct.gamma_g2 * self.sk // kind of a bls signature on gamma_g2
    }
}

impl AggregateKey {
    pub fn new(pk: Vec<PublicKey>, params: &PowersOfTau) -> Self {
        let n = pk.len();
        let h_minus1 = G2Projective::generator() * (-Scalar::ONE);
        let z_g2 = G2Projective::from(params.powers_of_h[n]) + h_minus1;

        // gather sk_li from all public keys
        let mut ask = G1Projective::identity();
        for pki in pk.iter() {
            ask += pki.sk_li;
        }

        let mut agg_sk_li_lj_z = vec![];
        for i in 0..n {
            let mut agg_sk_li_lj_zi = G1Projective::identity();
            for pkj in pk.iter() {
                agg_sk_li_lj_zi += pkj.sk_li_lj_z[i];
            }
            agg_sk_li_lj_z.push(agg_sk_li_lj_zi);
        }

        // Compute pairing e(g, h)
        use pairing::Engine;
        let e_gh = blstrs::Bls12::pairing(
            &G1Affine::from(params.powers_of_g[0]),
            &G2Affine::from(params.powers_of_h[0]),
        );

        AggregateKey {
            pk,
            agg_sk_li_lj_z,
            ask,
            z_g2,
            h_minus1,
            e_gh,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::OsRng;

    #[test]
    fn test_setup() {
        let mut rng = OsRng;
        let n = 16;
        let tau = Scalar::random(&mut rng);
        let params = KZG10::setup(n, tau).unwrap();
        let lagrange_params = LagrangePowers::new(tau, n);

        let mut sk: Vec<SecretKey> = Vec::new();
        let mut pk: Vec<PublicKey> = Vec::new();
        let mut lagrange_pk: Vec<PublicKey> = Vec::new();

        for i in 0..n {
            sk.push(SecretKey::new(&mut rng));
            pk.push(sk[i].get_pk(i, &params, n));
            lagrange_pk.push(sk[i].lagrange_get_pk(i, &lagrange_params, n));

            assert_eq!(pk[i].sk_li, lagrange_pk[i].sk_li);
            assert_eq!(pk[i].sk_li_minus0, lagrange_pk[i].sk_li_minus0);
            assert_eq!(pk[i].sk_li_x, lagrange_pk[i].sk_li_x);
            assert_eq!(pk[i].sk_li_lj_z, lagrange_pk[i].sk_li_lj_z);
        }

        let _ak = AggregateKey::new(pk, &params);
    }
}
