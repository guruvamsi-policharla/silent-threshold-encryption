use crate::{crs::CRS, setup::AggregateKey};
use ark_ec::{pairing::Pairing, PrimeGroup};
use ark_serialize::*;
use ark_std::UniformRand;
use std::ops::Mul;

use aes_gcm::{aead::Aead, Aes256Gcm, Key, KeyInit};
use hkdf::Hkdf;
use sha2::Sha256;

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct Ciphertext<E: Pairing> {
    pub gamma_g2: E::G2,
    pub sa1: [E::G1; 2],
    pub sa2: [E::G2; 6],
    // pub enc_key: PairingOutput<E>, // key to be used for encapsulation (linearly homomorphic)
    pub ct: Vec<u8>, //encrypted message
    pub t: usize,    //threshold
}

impl<E: Pairing> Ciphertext<E> {
    pub fn new(gamma_g2: E::G2, sa1: [E::G1; 2], sa2: [E::G2; 6], ct: Vec<u8>, t: usize) -> Self {
        Ciphertext {
            gamma_g2,
            sa1,
            sa2,
            ct,
            t,
        }
    }
}

/// t is the threshold for encryption and apk is the aggregated public key
pub fn encrypt<E: Pairing>(
    apk: &AggregateKey<E>,
    t: usize,
    crs: &CRS<E>,
    m: &[u8],
) -> Ciphertext<E> {
    let mut rng = ark_std::test_rng();
    let gamma = E::ScalarField::rand(&mut rng);
    let gamma_g2 = crs.powers_of_h[0] * gamma;

    let g = crs.powers_of_g[0];
    let h = crs.powers_of_h[0];

    let mut sa1 = [E::G1::generator(); 2];
    let mut sa2 = [E::G2::generator(); 6];

    let s = (0..5)
        .map(|_| E::ScalarField::rand(&mut rng))
        .collect::<Vec<_>>();

    // sa1[0] = s0*ask + s3*g^{tau^{t}} + s4*g
    sa1[0] = (apk.ask * s[0]) + (crs.powers_of_g[t] * s[3]) + (crs.powers_of_g[0] * s[4]);

    // sa1[1] = s2*g
    sa1[1] = g * s[2];

    // sa2[0] = s0*h + s2*gamma_g2
    sa2[0] = (h * s[0]) + (gamma_g2 * s[2]);

    // sa2[1] = s0*z_g2
    sa2[1] = apk.z_g2 * s[0];

    // sa2[2] = s0*h^tau + s1*h^{tau^2}
    sa2[2] = crs.powers_of_h[1] * s[0] + crs.powers_of_h[2] * s[1];

    // sa2[3] = s1*h
    sa2[3] = h * s[1];

    // sa2[4] = s3*h
    sa2[4] = h * s[3];

    // sa2[5] = s4*h^{tau}
    sa2[5] = (crs.powers_of_h[1]) * s[4];

    // enc_key = s4*e_gh
    let enc_key = apk.e_gh.mul(s[4]);
    let mut enc_key_bytes = Vec::new();
    enc_key.serialize_compressed(&mut enc_key_bytes).unwrap();

    // derive an encapsulation key from enc_key using an HKDF
    let hk = Hkdf::<Sha256>::new(None, &enc_key_bytes);
    let mut aes_key = [0u8; 32];
    let mut aes_nonce = [0u8; 12];
    hk.expand(&[1], &mut aes_key).unwrap();
    hk.expand(&[2], &mut aes_nonce).unwrap();

    // encrypt the message m using the derived key
    let aes_key: &Key<Aes256Gcm> = &aes_key.into();
    let cipher = Aes256Gcm::new(&aes_key);
    let ct = cipher.encrypt(&aes_nonce.into(), m).unwrap();

    Ciphertext {
        gamma_g2,
        sa1,
        sa2,
        ct,
        t,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        crs::CRS,
        setup::{LagPublicKey, SecretKey},
    };

    type E = ark_bls12_381::Bls12_381;
    type G1 = <E as Pairing>::G1;
    type G2 = <E as Pairing>::G2;

    #[test]
    fn test_encryption() {
        let mut rng = ark_std::test_rng();
        let n = 8;
        let crs = CRS::new(n, &mut rng);

        let msg = b"Hello, world!";

        let mut sk: Vec<SecretKey<E>> = Vec::new();
        let mut pk: Vec<LagPublicKey<E>> = Vec::new();

        for i in 0..n {
            sk.push(SecretKey::<E>::new(&mut rng));
            pk.push(sk[i].get_lagrange_pk(i, &crs))
        }

        let ak = AggregateKey::<E>::new(pk, &crs);
        let ct = encrypt::<E>(&ak, 2, &crs, msg);

        let mut ct_bytes = Vec::new();
        ct.serialize_compressed(&mut ct_bytes).unwrap();
        println!("Compressed ciphertext: {} bytes", ct_bytes.len());

        let mut g1_bytes = Vec::new();
        let mut g2_bytes = Vec::new();
        let mut e_gh_bytes = Vec::new();

        let g = G1::generator();
        let h = G2::generator();

        g.serialize_compressed(&mut g1_bytes).unwrap();
        h.serialize_compressed(&mut g2_bytes).unwrap();
        ak.e_gh.serialize_compressed(&mut e_gh_bytes).unwrap();

        println!("G1 len: {} bytes", g1_bytes.len());
        println!("G2 len: {} bytes", g2_bytes.len());
        println!("GT len: {} bytes", e_gh_bytes.len());
    }
}
