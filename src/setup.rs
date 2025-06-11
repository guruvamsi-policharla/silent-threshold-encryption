use crate::{
    crs::CRS,
    encryption::Ciphertext,
    utils::{lagrange_poly, open_all_values},
};
use ark_ec::{pairing::Pairing, AffineRepr, PrimeGroup, VariableBaseMSM};
use ark_ff::FftField;
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial,
    Radix2EvaluationDomain,
};
use ark_serialize::*;
use ark_std::{rand::RngCore, UniformRand, Zero};

use crate::utils::{ark_de, ark_se};
use serde::{Deserialize, Serialize};

#[derive(Clone)]
pub struct LagPolys<F: FftField> {
    pub l: Vec<DensePolynomial<F>>,
    pub l_minus0: Vec<DensePolynomial<F>>,
    pub l_x: Vec<DensePolynomial<F>>,
    pub li_lj_z: Vec<Vec<DensePolynomial<F>>>,
    pub denom: F,
}

impl<F: FftField> LagPolys<F> {
    // domain is the roots of unity of size n
    pub fn new(n: usize) -> Self {
        let domain = Radix2EvaluationDomain::<F>::new(n).unwrap();

        // compute polynomial L_i(X)
        let mut l = vec![DensePolynomial::zero(); n];
        for i in 0..n {
            l[i] = lagrange_poly(n, i);
        }

        // compute polynomial (L_i(X) - L_i(0))*X
        let mut l_minus0 = vec![DensePolynomial::zero(); n];
        for i in 0..n {
            let mut li_minus0_coeffs = l[i].coeffs.clone();
            li_minus0_coeffs[0] = F::zero();
            li_minus0_coeffs.insert(0, F::zero());
            l_minus0[i] = DensePolynomial::from_coefficients_vec(li_minus0_coeffs);
        }

        // compute polynomial (L_i(X) - L_i(0))/X
        let mut l_x = vec![DensePolynomial::zero(); n];
        for i in 0..n {
            l_x[i] = DensePolynomial::from_coefficients_vec(l_minus0[i].coeffs[2..].to_vec());
        }

        // compute polynomial L_i(X)*L_j(X)/Z(X) and (L_i(X)*L_i(X) - L_i(X))/Z(X)
        let mut li_lj_z = vec![vec![DensePolynomial::zero(); n]; n];
        for i in 0..n {
            for j in 0..n {
                li_lj_z[i][j] = if i == j {
                    (&l[i] * &l[i] - &l[i]).divide_by_vanishing_poly(domain).0
                } else {
                    (&l[i] * &l[j]).divide_by_vanishing_poly(domain).0
                };
            }
        }

        let mut denom = F::one();
        for i in 1..n {
            denom *= F::one() - domain.element(i);
        }

        // for i in 0..n {
        //     for j in 0..n {
        //         let monomial =
        //             DensePolynomial::from_coefficients_vec(vec![-domain.element(j), F::one()]);

        //         let computed = &l[i] / &monomial;
        //         assert_eq!(
        //             li_lj_z[i][j].evaluate(&F::zero()),
        //             computed.evaluate(&F::zero()) / (denom * domain.element(n - j))
        //         );
        //     }
        // }

        Self {
            l,
            l_minus0,
            l_x,
            li_lj_z,
            denom: denom.inverse().unwrap(),
        }
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Serialize, Deserialize, Clone)]
pub struct SecretKey<E: Pairing> {
    pub id: usize,
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    sk: E::ScalarField,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PartialDecryption<E: Pairing> {
    /// Party id
    pub id: usize,
    /// Party commitment
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub signature: E::G2,
}

impl<E: Pairing> PartialDecryption<E> {
    pub fn zero() -> Self {
        PartialDecryption {
            id: 0,
            signature: E::G2::zero(),
        }
    }
}

/// Position oblivious public key -- slower to aggregate
#[derive(CanonicalSerialize, CanonicalDeserialize, Serialize, Deserialize, Clone)]
pub struct PublicKey<E: Pairing> {
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub bls_pk: E::G1, //BLS pk
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub hints: Vec<E::G1Affine>, //hints
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub y: Vec<E::G1Affine>, /* preprocessed toeplitz matrix. only for efficiency and can be
                              * computed from hints */
    pub id: usize, // canonically assigned unique id in the system
}

/// Public key that can only be used in a fixed position -- faster to aggregate
#[derive(CanonicalSerialize, CanonicalDeserialize, Serialize, Deserialize, Clone)]
pub struct LagPublicKey<E: Pairing> {
    pub id: usize,       //id of the party
    pub position: usize, //position in the aggregate key
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub bls_pk: E::G1, //BLS pk
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub sk_li: E::G1, //hint
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub sk_li_minus0: E::G1, //hint
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub sk_li_lj_z: Vec<E::G1>, //hint
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub sk_li_x: E::G1, //hint
}

impl<E: Pairing> LagPublicKey<E> {
    pub fn new(
        id: usize,
        position: usize,
        bls_pk: E::G1,
        sk_li: E::G1,
        sk_li_minus0: E::G1,
        sk_li_lj_z: Vec<E::G1>, //i = id
        sk_li_x: E::G1,
    ) -> Self {
        LagPublicKey {
            id,
            position,
            bls_pk,
            sk_li,
            sk_li_minus0,
            sk_li_lj_z,
            sk_li_x,
        }
    }
}

impl<E: Pairing> SecretKey<E> {
    pub fn new<R: RngCore>(rng: &mut R, id: usize) -> Self {
        SecretKey {
            id,
            sk: E::ScalarField::rand(rng),
        }
    }

    pub fn from_scalar(sk: E::ScalarField, id: usize) -> Self {
        SecretKey { id, sk }
    }

    pub fn get_pk(&self, crs: &CRS<E>) -> PublicKey<E> {
        let mut hints = vec![E::G1Affine::zero(); crs.powers_of_g.len()];

        let bls_pk = E::G1::generator() * self.sk;

        for i in 0..crs.powers_of_g.len() {
            hints[i] = (crs.powers_of_g[i] * self.sk).into();
        }

        // compute y
        let mut y = vec![E::G1Affine::zero(); crs.y.len()];
        for i in 0..crs.y.len() {
            y[i] = (crs.y[i] * self.sk).into();
        }

        PublicKey {
            id: self.id,
            bls_pk,
            hints,
            y,
        }
    }

    pub fn get_lagrange_pk(&self, position: usize, crs: &CRS<E>) -> LagPublicKey<E> {
        let mut sk_li_lj_z = vec![];

        let sk_li = crs.li[position] * self.sk;

        let sk_li_minus0 = crs.li_minus0[position] * self.sk;

        let sk_li_x = crs.li_x[position] * self.sk;

        for j in 0..crs.n {
            sk_li_lj_z.push(crs.li_lj_z[position][j] * self.sk);
        }

        LagPublicKey {
            id: self.id,
            position,
            bls_pk: E::G1::generator() * self.sk,
            sk_li,
            sk_li_minus0,
            sk_li_lj_z,
            sk_li_x,
        }
    }

    pub fn partial_decryption(&self, ct: &Ciphertext<E>) -> PartialDecryption<E> {
        PartialDecryption {
            id: self.id,
            signature: ct.gamma_g2 * self.sk, // bls signature on gamma_g2
        }
    }
}

impl<E: Pairing> PublicKey<E> {
    pub fn get_lag_public_key(
        &self,
        position: usize,
        crs: &CRS<E>,
        lag_polys: &LagPolys<E::ScalarField>,
    ) -> LagPublicKey<E> {
        assert!(position < crs.n, "position out of bounds");

        let bls_pk = self.bls_pk;

        // compute sk_li
        let sk_li = E::G1::msm(
            &self.hints[0..lag_polys.l[position].degree() + 1],
            &lag_polys.l[position],
        )
        .unwrap();

        // compute sk_li_minus0
        let sk_li_minus0 = E::G1::msm(
            &self.hints[0..lag_polys.l_minus0[position].degree() + 1],
            &lag_polys.l_minus0[position],
        )
        .unwrap();

        // compute sk_li_x
        let sk_li_x = E::G1::msm(
            &self.hints[0..lag_polys.l_x[position].degree() + 1],
            &lag_polys.l_x[position],
        )
        .unwrap();

        // compute sk*Li*Lj/Z = sk*Li/(X-omega^j)*(omega^j/denom) for all j in [n]\{i}
        // for j = i: (Li^2 - Li)/Z = (Li - 1)/(X-omega^i)*(omega^i/denom)
        // this is the same as computing KZG opening proofs at all points
        // in the roots of unity domain for the polynomial Li(X), where the
        // crs is {g^sk, g^{sk * tau}, g^{sk * tau^2}, ...}
        // todo: move to https://eprint.iacr.org/2024/1279.pdf
        let domain = Radix2EvaluationDomain::<E::ScalarField>::new(crs.n).unwrap();
        let mut sk_li_lj_z = open_all_values::<E>(&self.y, &lag_polys.l[position].coeffs, &domain);
        for j in 0..crs.n {
            sk_li_lj_z[j] *= domain.element(j) * lag_polys.denom;
        }

        // // compute sk_li_lj_z
        // let mut sk_li_lj_z = vec![E::G1::zero(); crs.n];

        // let timer = start_timer!(|| "msm version");
        // for j in 0..crs.n {
        //     sk_li_lj_z[j] = E::G1::msm(
        //         &self.hints[0..lag_polys.li_lj_z[id][j].degree() + 1],
        //         &lag_polys.li_lj_z[id][j],
        //     )
        //     .unwrap();
        // }
        // end_timer!(timer);

        // assert_eq!(sk_li_lj_z, my_sk_li_lj_z);

        LagPublicKey {
            id: self.id,
            position,
            bls_pk,
            sk_li,
            sk_li_minus0,
            sk_li_lj_z,
            sk_li_x,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aggregate::AggregateKey;
    type E = ark_bls12_381::Bls12_381;
    type F = ark_bls12_381::Fr;

    #[test]
    fn test_setup() {
        let mut rng = ark_std::test_rng();
        let n = 1 << 4;
        let crs = CRS::<E>::new(n, &mut rng);

        let mut sk: Vec<SecretKey<E>> = Vec::new();
        let mut pk: Vec<LagPublicKey<E>> = Vec::new();
        let mut lagrange_pk: Vec<LagPublicKey<E>> = Vec::new();

        for i in 0..n {
            sk.push(SecretKey::<E>::new(&mut rng, i));
            pk.push(sk[i].get_lagrange_pk(i, &crs));
            lagrange_pk.push(sk[i].get_lagrange_pk(i, &crs));

            assert_eq!(pk[i].sk_li, lagrange_pk[i].sk_li);
            assert_eq!(pk[i].sk_li_minus0, lagrange_pk[i].sk_li_minus0);
            assert_eq!(pk[i].sk_li_x, lagrange_pk[i].sk_li_x); //computed incorrectly go fix it
            assert_eq!(pk[i].sk_li_lj_z, lagrange_pk[i].sk_li_lj_z);
        }

        let _ak = AggregateKey::<E>::new(pk, &crs);
    }

    #[test]
    fn test_setup_lag_setup() {
        let mut rng = ark_std::test_rng();
        let n = 1 << 7;
        let crs = CRS::<E>::new(n, &mut rng);
        let lagpolys = LagPolys::<F>::new(n);

        let sk = SecretKey::<E>::new(&mut rng, 0);
        let pk = sk.get_pk(&crs);
        let lag_pk = sk.get_lagrange_pk(0, &crs);

        let computed_lag_pk = pk.get_lag_public_key(0, &crs, &lagpolys);

        assert_eq!(computed_lag_pk.bls_pk, lag_pk.bls_pk);
        assert_eq!(computed_lag_pk.sk_li, lag_pk.sk_li);
        assert_eq!(computed_lag_pk.sk_li_minus0, lag_pk.sk_li_minus0);
        assert_eq!(computed_lag_pk.sk_li_x, lag_pk.sk_li_x);
        assert_eq!(computed_lag_pk.sk_li_lj_z, lag_pk.sk_li_lj_z);
    }
}
