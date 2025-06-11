use ark_ec::{pairing::Pairing, VariableBaseMSM};
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial,
    Radix2EvaluationDomain,
};
use ark_serialize::CanonicalSerialize;
use ark_std::{One, Zero};

use aes_gcm::{aead::Aead, Aes256Gcm, Key, KeyInit};
use hkdf::Hkdf;
use sha2::Sha256;

use crate::{
    aggregate::AggregateKey, crs::CRS, encryption::Ciphertext, setup::PartialDecryption,
    utils::interp_mostly_zero,
};

pub fn agg_dec<E: Pairing>(
    partial_decryptions: &Vec<PartialDecryption<E>>, /* insert 0 if a party did not respond or
                                                      * verification failed */
    ct: &Ciphertext<E>,
    selector: &[bool],
    agg_key: &AggregateKey<E>,
    crs: &CRS<E>,
) -> Vec<u8> {
    let domain = Radix2EvaluationDomain::<E::ScalarField>::new(crs.n).unwrap();
    let domain_elements: Vec<E::ScalarField> = domain.elements().collect();

    // points is where B is set to zero
    // parties is the set of parties who have signed
    let mut points = vec![];
    let mut parties: Vec<usize> = Vec::new(); // parties indexed from 0..n-1
    for i in 0..crs.n {
        if selector[i] {
            parties.push(i);
        } else {
            points.push(domain_elements[i]);
        }
    }

    let b = interp_mostly_zero(&points);
    let b_evals = domain.fft(&b.coeffs);

    debug_assert_eq!(
        b.degree(),
        points.len(),
        "b.degree should be equal to points.len()"
    );
    debug_assert!(b.evaluate(&E::ScalarField::zero()) == E::ScalarField::one());

    // commit to b in g2
    let b_g2: E::G2 = crs.commit_g2(&b.coeffs);

    // q0 = (b-1)/x
    let q0_g1 = crs.compute_opening_proof(&b.coeffs, &E::ScalarField::zero());

    // bhat = x^{t} * b
    // insert t 0s at the beginning of bhat.coeffs
    let mut bhat_coeffs = vec![E::ScalarField::zero(); ct.t];
    bhat_coeffs.append(&mut b.coeffs.clone());
    let bhat = DensePolynomial::from_coefficients_vec(bhat_coeffs);
    debug_assert_eq!(bhat.degree(), crs.n);

    let bhat_g1: E::G1 = crs.commit_g1(&bhat.coeffs);

    let n_inv = E::ScalarField::one() / E::ScalarField::from(crs.n as u32);

    // compute the aggregate public key
    let mut bases: Vec<<E as Pairing>::G1Affine> = Vec::new();
    let mut scalars: Vec<<E as Pairing>::ScalarField> = Vec::new();
    for &i in &parties {
        bases.push(agg_key.lag_pks[i].bls_pk.into());
        scalars.push(b_evals[i]);
    }
    let mut apk = E::G1::msm(bases.as_slice(), scalars.as_slice()).unwrap();
    apk *= n_inv;

    // compute sigma = (\sum B(omega^i)partial_decryptions[i])/(n) for i in parties
    let mut bases: Vec<<E as Pairing>::G2Affine> = Vec::new();
    let mut scalars: Vec<<E as Pairing>::ScalarField> = Vec::new();
    for &i in &parties {
        bases.push(partial_decryptions[i].signature.into());
        scalars.push(b_evals[i]);
    }
    let mut sigma = E::G2::msm(bases.as_slice(), scalars.as_slice()).unwrap();
    sigma *= n_inv;

    // compute Qx, Qhatx and Qz
    let mut bases: Vec<<E as Pairing>::G1Affine> = Vec::new();
    let mut scalars: Vec<<E as Pairing>::ScalarField> = Vec::new();
    for &i in &parties {
        bases.push(agg_key.lag_pks[i].sk_li_x.into());
        scalars.push(b_evals[i]);
    }
    let qx = E::G1::msm(bases.as_slice(), scalars.as_slice()).unwrap();

    let mut bases: Vec<<E as Pairing>::G1Affine> = Vec::new();
    let mut scalars: Vec<<E as Pairing>::ScalarField> = Vec::new();
    for &i in &parties {
        bases.push(agg_key.lag_pks[i].sk_li_minus0.into());
        scalars.push(b_evals[i]);
    }
    let qhatx = E::G1::msm(bases.as_slice(), scalars.as_slice()).unwrap();

    let mut bases: Vec<<E as Pairing>::G1Affine> = Vec::new();
    let mut scalars: Vec<<E as Pairing>::ScalarField> = Vec::new();
    for &i in &parties {
        bases.push(agg_key.agg_sk_li_lj_z[i].into());
        scalars.push(b_evals[i]);
    }
    let qz = E::G1::msm(bases.as_slice(), scalars.as_slice()).unwrap();

    // e(w1||sa1, sa2||w2)
    let minus1 = -E::ScalarField::one();
    let w1 = [
        apk * (minus1),
        qz * (minus1),
        qx * (minus1),
        qhatx,
        bhat_g1 * (minus1),
        q0_g1 * (minus1),
    ];
    let w2 = [b_g2, sigma];

    let mut enc_key_lhs = w1.to_vec();
    enc_key_lhs.append(&mut ct.sa1.to_vec());

    let mut enc_key_rhs = ct.sa2.to_vec();
    enc_key_rhs.append(&mut w2.to_vec());

    let enc_key = E::multi_pairing(enc_key_lhs, enc_key_rhs);
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

    cipher.decrypt(&aes_nonce.into(), ct.ct.as_ref()).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        crs::CRS,
        encryption::encrypt,
        setup::{PartialDecryption, SecretKey},
    };

    type E = ark_bls12_381::Bls12_381;
    type G2 = <E as Pairing>::G2;
    use ark_std::UniformRand;

    #[test]
    fn test_decryption() {
        let mut rng = ark_std::test_rng();
        let n = 1 << 3; // actually n-1 total parties. one party is a dummy party that is always true
        let t: usize = 1;
        debug_assert!(t < n);

        let crs = CRS::new(n, &mut rng);

        let msg = b"Hello, world!";

        let sk = (0..n)
            .map(|i| SecretKey::<E>::new(&mut rng, i))
            .collect::<Vec<_>>();

        let pk = sk
            .iter()
            .enumerate()
            .map(|(i, sk)| sk.get_lagrange_pk(i, &crs))
            .collect::<Vec<_>>();

        let (ak, ek) = AggregateKey::<E>::new(pk, &crs);

        let gamma_g2 = G2::rand(&mut rng);
        let ct = encrypt::<E>(&ek, t, &crs, gamma_g2, msg);

        // compute partial decryptions
        let mut partial_decryptions: Vec<PartialDecryption<E>> = Vec::new();
        for i in 0..t {
            partial_decryptions.push(sk[i].partial_decryption(&ct));
        }
        for _ in t..n {
            partial_decryptions.push(PartialDecryption::<E>::zero());
        }

        // compute the decryption key
        let mut selector: Vec<bool> = Vec::new();
        for _ in 0..t {
            selector.push(true);
        }
        for _ in t..n {
            selector.push(false);
        }

        assert_eq!(
            agg_dec(&partial_decryptions, &ct, &selector, &ak, &crs),
            msg
        );
    }
}
