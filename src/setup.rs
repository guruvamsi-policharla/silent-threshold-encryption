use crate::crs::CRS;
use crate::encryption::Ciphertext;
use crate::utils::lagrange_poly;
use ark_ec::{
    pairing::{Pairing, PairingOutput},
    AffineRepr, PrimeGroup, VariableBaseMSM,
};
use ark_ff::FftField;
use ark_poly::DenseUVPolynomial;
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, Polynomial, Radix2EvaluationDomain};
use ark_serialize::*;
use ark_std::{rand::RngCore, One, UniformRand, Zero};

#[derive(Clone)]
pub struct LagPolys<F: FftField> {
    pub l: Vec<DensePolynomial<F>>,
    pub l_minus0: Vec<DensePolynomial<F>>,
    pub l_x: Vec<DensePolynomial<F>>,
    pub li_lj_z: Vec<Vec<DensePolynomial<F>>>,
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

        Self {
            l,
            l_minus0,
            l_x,
            li_lj_z,
        }
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct SecretKey<E: Pairing> {
    sk: E::ScalarField,
}

/// Position oblivious public key -- slower to aggregate
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct PublicKey<E: Pairing> {
    pub bls_pk: E::G1,           //BLS pk
    pub hints: Vec<E::G1Affine>, //hints
}

/// Public key that can only be used in a fixed position -- faster to aggregate
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct LagPublicKey<E: Pairing> {
    pub id: usize,
    pub bls_pk: E::G1,          //BLS pk
    pub sk_li: E::G1,           //hint
    pub sk_li_minus0: E::G1,    //hint
    pub sk_li_lj_z: Vec<E::G1>, //hint
    pub sk_li_x: E::G1,         //hint
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct AggregateKey<E: Pairing> {
    pub pk: Vec<LagPublicKey<E>>,
    pub agg_sk_li_lj_z: Vec<E::G1>,
    pub ask: E::G1,
    pub z_g2: E::G2,

    //preprocessed values
    pub e_gh: PairingOutput<E>,
}

impl<E: Pairing> LagPublicKey<E> {
    pub fn new(
        id: usize,
        bls_pk: E::G1,
        sk_li: E::G1,
        sk_li_minus0: E::G1,
        sk_li_lj_z: Vec<E::G1>, //i = id
        sk_li_x: E::G1,
    ) -> Self {
        LagPublicKey {
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

    pub fn from_scalar(sk: E::ScalarField) -> Self {
        SecretKey { sk }
    }

    pub fn get_pk(&self, crs: &CRS<E>) -> PublicKey<E> {
        let mut hints = vec![E::G1Affine::zero(); crs.powers_of_g.len()];

        let bls_pk = E::G1::generator() * self.sk;

        for i in 0..crs.powers_of_g.len() {
            hints[i] = (crs.powers_of_g[i] * self.sk).into();
        }

        PublicKey { bls_pk, hints }
    }

    pub fn get_lagrange_pk(&self, id: usize, crs: &CRS<E>) -> LagPublicKey<E> {
        let mut sk_li_lj_z = vec![];

        let sk_li = crs.li[id] * self.sk;

        let sk_li_minus0 = crs.li_minus0[id] * self.sk;

        let sk_li_x = crs.li_x[id] * self.sk;

        for j in 0..crs.n {
            sk_li_lj_z.push(crs.li_lj_z[id][j] * self.sk);
        }

        LagPublicKey {
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

impl<E: Pairing> PublicKey<E> {
    pub fn get_lag_public_key(
        &self,
        id: usize,
        pk: &PublicKey<E>,
        crs: &CRS<E>,
        lagpolys: &LagPolys<E::ScalarField>,
    ) -> LagPublicKey<E> {
        let bls_pk = pk.bls_pk;

        // compute sk_li
        let sk_li = E::G1::msm(&pk.hints[0..lagpolys.l[id].degree() + 1], &lagpolys.l[id]).unwrap();

        // compute sk_li_minus0
        let sk_li_minus0 = E::G1::msm(
            &pk.hints[0..lagpolys.l_minus0[id].degree() + 1],
            &lagpolys.l_minus0[id],
        )
        .unwrap();

        // compute sk_li_x
        let sk_li_x = E::G1::msm(
            &pk.hints[0..lagpolys.l_x[id].degree() + 1],
            &lagpolys.l_x[id],
        )
        .unwrap();

        // compute sk_li_lj_z
        let mut sk_li_lj_z = vec![E::G1::zero(); crs.n];

        for j in 0..crs.n {
            sk_li_lj_z[j] = E::G1::msm(
                &pk.hints[0..lagpolys.li_lj_z[id][j].degree() + 1],
                &lagpolys.li_lj_z[id][j],
            )
            .unwrap();
        }

        LagPublicKey {
            id,
            bls_pk,
            sk_li,
            sk_li_minus0,
            sk_li_lj_z,
            sk_li_x,
        }
    }
}

impl<E: Pairing> AggregateKey<E> {
    pub fn new(pk: Vec<LagPublicKey<E>>, crs: &CRS<E>) -> Self {
        let n = pk.len();
        let z_g2 = crs.powers_of_h[n] + crs.powers_of_h[0] * (-E::ScalarField::one());

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
            e_gh: E::pairing(crs.powers_of_g[0], crs.powers_of_h[0]),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
            sk.push(SecretKey::<E>::new(&mut rng));
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
        let n = 1 << 4;
        let crs = CRS::<E>::new(n, &mut rng);
        let lagpolys = LagPolys::<F>::new(n);

        let sk = SecretKey::<E>::new(&mut rng);
        let pk = sk.get_pk(&crs);
        let lag_pk = sk.get_lagrange_pk(0, &crs);

        let computed_lag_pk = pk.get_lag_public_key(0, &pk, &crs, &lagpolys);
        assert_eq!(computed_lag_pk.bls_pk, lag_pk.bls_pk);
        assert_eq!(computed_lag_pk.sk_li, lag_pk.sk_li);
        assert_eq!(computed_lag_pk.sk_li_minus0, lag_pk.sk_li_minus0);
        assert_eq!(computed_lag_pk.sk_li_x, lag_pk.sk_li_x);
        assert_eq!(computed_lag_pk.sk_li_lj_z, lag_pk.sk_li_lj_z);
    }
}
