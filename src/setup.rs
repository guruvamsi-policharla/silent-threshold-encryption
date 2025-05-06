use crate::crs::CRS;
use crate::encryption::Ciphertext;
use ark_ec::pairing::PairingOutput;
use ark_ec::{pairing::Pairing, PrimeGroup};
use ark_serialize::*;
use ark_std::{rand::RngCore, One, UniformRand, Zero};

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

    pub fn get_pk(&self, id: usize, crs: &CRS<E>, n: usize) -> PublicKey<E> {
        let mut sk_li_lj_z = vec![];

        let sk_li = crs.li[id] * self.sk;

        let sk_li_minus0 = crs.li_minus0[id] * self.sk;

        let sk_li_x = crs.li_x[id] * self.sk;

        for j in 0..n {
            sk_li_lj_z.push(crs.li_lj_z[id][j] * self.sk);
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
    pub fn new(pk: Vec<PublicKey<E>>, crs: &CRS<E>) -> Self {
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

    #[test]
    fn test_setup() {
        let mut rng = ark_std::test_rng();
        let n = 16;
        let crs = CRS::<E>::new(n, &mut rng);

        let mut sk: Vec<SecretKey<E>> = Vec::new();
        let mut pk: Vec<PublicKey<E>> = Vec::new();
        let mut lagrange_pk: Vec<PublicKey<E>> = Vec::new();

        for i in 0..n {
            sk.push(SecretKey::<E>::new(&mut rng));
            pk.push(sk[i].get_pk(i, &crs, n));
            lagrange_pk.push(sk[i].get_pk(i, &crs, n));

            assert_eq!(pk[i].sk_li, lagrange_pk[i].sk_li);
            assert_eq!(pk[i].sk_li_minus0, lagrange_pk[i].sk_li_minus0);
            assert_eq!(pk[i].sk_li_x, lagrange_pk[i].sk_li_x); //computed incorrectly go fix it
            assert_eq!(pk[i].sk_li_lj_z, lagrange_pk[i].sk_li_lj_z);
        }

        let _ak = AggregateKey::<E>::new(pk, &crs);
    }
}
