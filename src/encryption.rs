use crate::{kzg::PowersOfTau, setup::AggregateKey};
use blstrs::{G1Projective, G2Projective, Gt, Scalar};
use ff::Field;
use group::Group;
use serde::{Deserialize, Serialize};

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct Ciphertext {
    pub gamma_g2: G2Projective,
    pub sa1: [G1Projective; 2],
    pub sa2: [G2Projective; 6],
    pub enc_key: Gt, //key to be used for encapsulation
    pub t: usize,    //threshold
}

impl Ciphertext {
    pub fn new(
        gamma_g2: G2Projective,
        sa1: [G1Projective; 2],
        sa2: [G2Projective; 6],
        enc_key: Gt,
        t: usize,
    ) -> Self {
        Ciphertext {
            gamma_g2,
            sa1,
            sa2,
            enc_key,
            t,
        }
    }
}

/// t is the threshold for encryption and apk is the aggregated public key
pub fn encrypt(apk: &AggregateKey, t: usize, params: &PowersOfTau) -> Ciphertext {
    use rand_core::OsRng;
    let mut rng = OsRng;

    let gamma = Scalar::random(&mut rng);
    let gamma_g2 = G2Projective::from(params.powers_of_h[0]) * gamma;

    let g = G1Projective::from(params.powers_of_g[0]);
    let h = G2Projective::from(params.powers_of_h[0]);

    let mut sa1 = [G1Projective::generator(); 2];
    let mut sa2 = [G2Projective::generator(); 6];

    let mut s: [Scalar; 5] = [Scalar::ZERO; 5];

    s.iter_mut().for_each(|s| *s = Scalar::random(&mut rng));

    // sa1[0] = s0*ask + s3*g^{tau^{t+1}} + s4*g
    sa1[0] = (apk.ask * s[0])
        + (G1Projective::from(params.powers_of_g[t + 1]) * s[3])
        + (G1Projective::from(params.powers_of_g[0]) * s[4]);

    // sa1[1] = s2*g
    sa1[1] = g * s[2];

    // sa2[0] = s0*h + s2*gamma_g2
    sa2[0] = (h * s[0]) + (gamma_g2 * s[2]);

    // sa2[1] = s0*z_g2
    sa2[1] = apk.z_g2 * s[0];

    // sa2[2] = s0*h^tau + s1*h^tau
    sa2[2] = G2Projective::from(params.powers_of_h[1]) * (s[0] + s[1]);

    // sa2[3] = s1*h
    sa2[3] = h * s[1];

    // sa2[4] = s3*h
    sa2[4] = h * s[3];

    // sa2[5] = s4*h^{tau - omega^0}
    sa2[5] = (G2Projective::from(params.powers_of_h[1]) + apk.h_minus1) * s[4];

    // enc_key = e_gh^s4
    let enc_key = apk.e_gh * s[4];

    Ciphertext {
        gamma_g2,
        sa1,
        sa2,
        enc_key,
        t,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        kzg::KZG10,
        setup::{PublicKey, SecretKey},
    };
    use rand_core::OsRng;

    #[test]
    fn test_encryption() {
        let mut rng = OsRng;
        let n = 8;
        let tau = Scalar::random(&mut rng);
        let params = KZG10::setup(n, tau).unwrap();

        let mut sk: Vec<SecretKey> = Vec::new();
        let mut pk: Vec<PublicKey> = Vec::new();

        for i in 0..n {
            sk.push(SecretKey::new(&mut rng));
            pk.push(sk[i].get_pk(0, &params, n))
        }

        let ak = AggregateKey::new(pk, &params);
        let ct = encrypt(&ak, 2, &params);

        // Test serialization sizes
        let ct_bytes = bincode::serialize(&ct).unwrap();
        println!("Ciphertext: {} bytes", ct_bytes.len());

        let g = G1Projective::generator();
        let h = G2Projective::generator();

        let g_bytes = bincode::serialize(&g).unwrap();
        let h_bytes = bincode::serialize(&h).unwrap();
        let e_gh_bytes = bincode::serialize(&ak.e_gh).unwrap();

        println!("G1 len: {} bytes", g_bytes.len());
        println!("G2 len: {} bytes", h_bytes.len());
        println!("GT len: {} bytes", e_gh_bytes.len());
    }
}
