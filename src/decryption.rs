use ark_ec::{
    pairing::{Pairing, PairingOutput},
    VariableBaseMSM,
};
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial,
    Radix2EvaluationDomain,
};
use ark_std::{One, Zero};
use std::ops::Div;

use crate::{
    encryption::Ciphertext,
    kzg::{PowersOfTau, KZG10},
    setup::AggregateKey,
    utils::interp_mostly_zero,
};

pub fn agg_dec<E: Pairing>(
    partial_decryptions: &[E::G2], //insert 0 if a party did not respond or verification failed
    ct: &Ciphertext<E>,
    selector: &[bool],
    agg_key: &AggregateKey<E>,
    params: &PowersOfTau<E>,
) -> PairingOutput<E> {
    let n = agg_key.pk.len();
    let domain = Radix2EvaluationDomain::<E::ScalarField>::new(n).unwrap();
    let domain_elements: Vec<E::ScalarField> = domain.elements().collect();

    // points is where B is set to zero
    // parties is the set of parties who have signed
    let mut points = vec![domain_elements[0]]; // 0 is the dummy party that is always true
    let mut parties: Vec<usize> = Vec::new(); // parties indexed from 0..n-1
    for i in 0..n {
        if selector[i] {
            parties.push(i);
        } else {
            points.push(domain_elements[i]);
        }
    }

    let b = interp_mostly_zero(E::ScalarField::one(), &points);
    let b_evals = domain.fft(&b.coeffs);

    debug_assert!(b.degree() == points.len() - 1);
    debug_assert!(b.evaluate(&domain_elements[0]) == E::ScalarField::one());

    // commit to b in g2
    let b_g2: E::G2 = KZG10::<E, DensePolynomial<E::ScalarField>>::commit_g2(params, &b)
        .unwrap()
        .into();

    // q0 = (b-1)/(x-domain_elements[0])
    let mut bminus1 = b.clone();
    bminus1.coeffs[0] -= E::ScalarField::one();

    debug_assert!(bminus1.evaluate(&domain_elements[0]) == E::ScalarField::zero());

    let xminus1 =
        DensePolynomial::from_coefficients_vec(vec![-domain_elements[0], E::ScalarField::one()]);
    let q0 = bminus1.div(&xminus1);

    let q0_g1: E::G1 = KZG10::<E, DensePolynomial<E::ScalarField>>::commit_g1(params, &q0)
        .unwrap()
        .into();

    // bhat = x^{t+1} * b
    // insert t+1 0s at the beginning of bhat.coeffs
    let mut bhat_coeffs = vec![E::ScalarField::zero(); ct.t + 1];
    bhat_coeffs.append(&mut b.coeffs.clone());
    let bhat = DensePolynomial::from_coefficients_vec(bhat_coeffs);
    debug_assert_eq!(bhat.degree(), n);

    let bhat_g1: E::G1 = KZG10::<E, DensePolynomial<E::ScalarField>>::commit_g1(params, &bhat)
        .unwrap()
        .into();

    let n_inv = E::ScalarField::one() / E::ScalarField::from((n) as u32);

    // compute the aggregate public key
    let mut bases: Vec<<E as Pairing>::G1Affine> = Vec::new();
    let mut scalars: Vec<<E as Pairing>::ScalarField> = Vec::new();
    for &i in &parties {
        bases.push(agg_key.pk[i].bls_pk.into());
        scalars.push(b_evals[i]);
    }
    let mut apk = E::G1::msm(bases.as_slice(), scalars.as_slice()).unwrap();
    apk *= n_inv;

    // compute sigma = (\sum B(omega^i)partial_decryptions[i])/(n) for i in parties
    let mut bases: Vec<<E as Pairing>::G2Affine> = Vec::new();
    let mut scalars: Vec<<E as Pairing>::ScalarField> = Vec::new();
    for &i in &parties {
        bases.push(partial_decryptions[i].into());
        scalars.push(b_evals[i]);
    }
    let mut sigma = E::G2::msm(bases.as_slice(), scalars.as_slice()).unwrap();
    sigma *= n_inv;

    // compute Qx, Qhatx and Qz
    let mut bases: Vec<<E as Pairing>::G1Affine> = Vec::new();
    let mut scalars: Vec<<E as Pairing>::ScalarField> = Vec::new();
    for &i in &parties {
        bases.push(agg_key.pk[i].sk_li_x.into());
        scalars.push(b_evals[i]);
    }
    let qx = E::G1::msm(bases.as_slice(), scalars.as_slice()).unwrap();

    let mut bases: Vec<<E as Pairing>::G1Affine> = Vec::new();
    let mut scalars: Vec<<E as Pairing>::ScalarField> = Vec::new();
    for &i in &parties {
        bases.push(agg_key.agg_sk_li_lj_z[i].into());
        scalars.push(b_evals[i]);
    }
    let qz = E::G1::msm(bases.as_slice(), scalars.as_slice()).unwrap();

    let mut bases: Vec<<E as Pairing>::G1Affine> = Vec::new();
    let mut scalars: Vec<<E as Pairing>::ScalarField> = Vec::new();
    for &i in &parties {
        bases.push(agg_key.pk[i].sk_li_minus0.into());
        scalars.push(b_evals[i]);
    }
    let qhatx = E::G1::msm(bases.as_slice(), scalars.as_slice()).unwrap();

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

    assert_eq!(enc_key, ct.enc_key);

    enc_key
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        encryption::encrypt,
        kzg::KZG10,
        setup::{PublicKey, SecretKey},
    };
    use ark_poly::univariate::DensePolynomial;
    use ark_std::UniformRand;

    type E = ark_bls12_381::Bls12_381;
    type G2 = <E as Pairing>::G2;
    type Fr = <E as Pairing>::ScalarField;
    type UniPoly381 = DensePolynomial<<E as Pairing>::ScalarField>;

    #[test]
    fn test_decryption() {
        let mut rng = ark_std::test_rng();
        let n = 1 << 4; // actually n-1 total parties. one party is a dummy party that is always true
        let t: usize = n / 2;
        debug_assert!(t < n);

        let tau = Fr::rand(&mut rng);
        let params = KZG10::<E, UniPoly381>::setup(n, tau.clone()).unwrap();

        let mut sk: Vec<SecretKey<E>> = Vec::new();
        let mut pk: Vec<PublicKey<E>> = Vec::new();

        // create the dummy party's keys
        sk.push(SecretKey::<E>::new(&mut rng));
        sk[0].nullify();
        pk.push(sk[0].get_pk(0, &params, n));

        for i in 1..n {
            sk.push(SecretKey::<E>::new(&mut rng));
            pk.push(sk[i].get_pk(i, &params, n))
        }

        let agg_key = AggregateKey::<E>::new(pk, &params);
        let ct = encrypt::<E>(&agg_key, t, &params);

        // compute partial decryptions
        let mut partial_decryptions: Vec<G2> = Vec::new();
        for i in 0..t + 1 {
            partial_decryptions.push(sk[i].partial_decryption(&ct));
        }
        for _ in t + 1..n {
            partial_decryptions.push(G2::zero());
        }

        // compute the decryption key
        let mut selector: Vec<bool> = Vec::new();
        for _ in 0..t + 1 {
            selector.push(true);
        }
        for _ in t + 1..n {
            selector.push(false);
        }

        let _dec_key = agg_dec(&partial_decryptions, &ct, &selector, &agg_key, &params);
    }
}
