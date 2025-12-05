use blstrs::{G1Projective, G2Projective, Gt, Scalar};
use ff::Field;
use group::Group;
use serde::{Deserialize, Serialize};

use crate::encryption::Ciphertext;
use crate::kzg::{PowersOfTau, KZG10};
use crate::polynomial::{DensePolynomial, Radix2EvaluationDomain};
use crate::utils::{lagrange_poly, lagrange_polys};
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
        let lagranges = lagrange_polys(n);
        let mut li_evals: Vec<Scalar> = vec![Scalar::ZERO; n];
        let mut li_evals_minus0: Vec<Scalar> = vec![Scalar::ZERO; n];
        let mut li_evals_x: Vec<Scalar> = vec![Scalar::ZERO; n];
        let tau_inv = tau.invert().unwrap();
        for (i, li) in lagranges.iter().enumerate() {
            li_evals[i] = li.evaluate(&tau);

            li_evals_minus0[i] = li_evals[i] - li.coeffs[0];

            li_evals_x[i] = li_evals_minus0[i] * tau_inv;
        }

        let z_eval = tau.pow([n as u64, 0, 0, 0]) - Scalar::ONE;
        let z_eval_inv = z_eval.invert().unwrap();

        let li = (0..n)
            .into_par_iter()
            .map(|i| G1Projective::generator() * li_evals[i])
            .collect::<Vec<_>>();

        let li_minus0 = li_evals_minus0
            .par_iter()
            .map(|&eval| G1Projective::generator() * eval)
            .collect::<Vec<_>>();

        let li_x = li_evals_x
            .par_iter()
            .map(|&eval| G1Projective::generator() * eval)
            .collect::<Vec<_>>();

        let li_lj_z = (0..n)
            .into_par_iter()
            .map(|i| {
                (0..n)
                    .into_par_iter()
                    .map(|j| {
                        if i == j {
                            G1Projective::generator()
                                * ((li_evals[i] * li_evals[i] - li_evals[i]) * z_eval_inv)
                        } else {
                            G1Projective::generator() * (li_evals[i] * li_evals[j] * z_eval_inv)
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

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
        let sk_li = params.li[id] * self.sk;

        let sk_li_minus0 = params.li_minus0[id] * self.sk;

        let sk_li_x = params.li_x[id] * self.sk;

        let sk_li_lj_z = params.li_lj_z[id]
            .iter()
            .take(n)
            .map(|elem| elem * self.sk)
            .collect();

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
        let h_minus1 = -G2Projective::generator();
        let z_g2 = G2Projective::from(params.powers_of_h[n]) + h_minus1;

        // gather sk_li from all public keys
        let mut ask = G1Projective::identity();
        for pki in pk.iter() {
            ask += pki.sk_li;
        }

        let agg_sk_li_lj_z = (0..n)
            .into_par_iter()
            .map(|i| {
                let mut agg_sk_li_lj_zi = G1Projective::identity();
                for pkj in pk.iter() {
                    agg_sk_li_lj_zi += pkj.sk_li_lj_z[i];
                }
                agg_sk_li_lj_zi
            })
            .collect::<Vec<_>>();

        AggregateKey {
            pk,
            agg_sk_li_lj_z,
            ask,
            z_g2,
            h_minus1,
            e_gh: params.e_gh,
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
