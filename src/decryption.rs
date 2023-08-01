use std::ops::Div;

use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ff::Field;
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain, Polynomial,
};
use ark_std::{One, Zero};

use crate::{
    encryption::Ciphertext,
    kzg::{UniversalParams, KZG10},
    setup::AggregateKey,
    utils::interp_mostly_zero,
};

pub fn agg_dec<E: Pairing>(
    partial_decryptions: &Vec<E::G2>, //insert 0 if a party did not respond or verification failed
    ct: &Ciphertext<E>,
    selector: &Vec<bool>,
    agg_key: &AggregateKey<E>,
    params: &UniversalParams<E>,
) -> PairingOutput<E> {
    let n = agg_key.pk.len();
    let domain = Radix2EvaluationDomain::<E::ScalarField>::new(n).unwrap();
    let domain_elements: Vec<E::ScalarField> = domain.elements().collect();

    // select the points that have selector[i] = false
    let mut points: Vec<E::ScalarField> = vec![domain_elements[0]]; // 0 is the dummy party that is always true
    let mut parties: Vec<usize> = Vec::new(); // parties indexed from 0..n-1
    selector.iter().enumerate().for_each(|(i, x)| {
        if *x==false {
            // println!("i:{}", i);
            points.push(domain_elements[i]);
            parties.push(i);
        }
    });

    let b = interp_mostly_zero(E::ScalarField::one(), &points);
    let b_evals = domain.fft(&b.coeffs);

    debug_assert!(b.degree() == points.len()-1);
    // println!("t:{}, points.len():{}", ct.t, points.len());

    // commit to b in g2
    let b_g2 = KZG10::<E, DensePolynomial<E::ScalarField>>::commit_g2(&params, &b)
        .unwrap()
        .into();

    // bhat = x^t * b
    // insert t 0s at the beginning of bhat.coeffs
    let mut bhat_coeffs = vec![E::ScalarField::zero(); ct.t];
    bhat_coeffs.append(&mut b.coeffs.clone());
    let bhat = DensePolynomial::from_coefficients_vec(bhat_coeffs);
    debug_assert_eq!(bhat.degree(),n-1);
    // println!("bhat.degree():{}", bhat.degree());

    let bhat = KZG10::<E, DensePolynomial<E::ScalarField>>::commit_g1(
        &params,
        &bhat,
    )
    .unwrap()
    .into();

    // todo: stop accumulating the dummy party public key
    // accumulate the public key
    // count number of true in selector
    let mut v: usize = 0;
    for i in 1..selector.len() {
        if selector[i] {
            v += 1;
        }
    }
    // selector.iter().for_each(|x| {
    //     if *x {
    //         v += 1;
    //     }
    // });
    let mut v = E::ScalarField::from(v as u32);
    v.inverse_in_place();

    // compute the aggregate public key
    let mut apk = E::G1::zero();
    for i in 1..agg_key.pk.len() {
        if selector[i] {
            apk += agg_key.pk[i].bls_pk;
        }
    }
    apk *= v;

    // compute Qx, Qhatx and Qz
    let mut qx: E::G1 = E::G1::zero();
    for &i in &parties {
        qx += agg_key.pk[i].sk_li_by_tau * b_evals[i];
    }

    let mut qz: E::G1 = E::G1::zero();
    for &i in &parties {
        let mut qz_i = E::G1::zero();
        for &j in &parties {
            qz_i += agg_key.pk[i].sk_li_by_z[j];
        }
        qz += qz_i * b_evals[i];
    }

    let mut qhatx: E::G1 = E::G1::zero();
    for &i in &parties {
        qhatx += agg_key.pk[i].sk_li_minus0 * b_evals[i];
    }

    // compute sigma = (\sum partial_decryptions[i])/v for i in parties
    let mut sigma: E::G2 = E::G2::zero();
    for &i in &parties {
        sigma += partial_decryptions[i];
    }
    sigma *= v;

    // q0 = (b-1)/(x-domain_elements[0])
    let mut bminus1 = b.clone();
    bminus1.coeffs[0] -= E::ScalarField::one();

    let q0 = b.div(&DensePolynomial::from_coefficients_vec(vec![
        -domain_elements[0],
        E::ScalarField::one(),
    ]));

    // debug_assert_eq!()

    let q0_g1 = KZG10::<E, DensePolynomial<E::ScalarField>>::commit_g1(&params, &q0)
        .unwrap()
        .into();

    let w1 = [apk, qz, qx, qhatx, bhat, q0_g1];
    let w2 = [sigma, b_g2];

    let mut enc_key_lhs = w1.to_vec();
    enc_key_lhs.append(&mut ct.sa1.to_vec());

    let mut enc_key_rhs = w2.to_vec();
    enc_key_rhs.append(&mut ct.sa2.to_vec());

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

    type E = ark_bls12_381::Bls12_381;
    // type F = <E as Pairing>::ScalarField;
    // type G1 = <E as Pairing>::G1;
    type G2 = <E as Pairing>::G2;
    type UniPoly381 = DensePolynomial<<E as Pairing>::ScalarField>;

    #[test]
    fn test_decryption() {
        let mut rng = ark_std::test_rng();
        let n = 1<<3; // actually n-1 total parties. one party is a dummy party that is always true
        let t: usize = 7;
        debug_assert!(t<n);

        let params = KZG10::<E, UniPoly381>::setup(n, &mut rng).unwrap();

        let mut sk: Vec<SecretKey<E>> = Vec::new();
        let mut pk: Vec<PublicKey<E>> = Vec::new();

        for i in 0..n {
            sk.push(SecretKey::<E>::new(&mut rng));
            pk.push(sk[i].get_pk(0, &params, n))
        }

        let agg_key = AggregateKey::<E>::new(pk, &params);
        let ct = encrypt::<E>(&agg_key, t, &params);

        // compute partial decryptions
        let mut partial_decryptions: Vec<G2> = Vec::new();
        for i in 0..t+1 {
            partial_decryptions.push(sk[i].partial_decryption(&ct));
        }
        for _ in t+1..n {
            partial_decryptions.push(G2::zero());
        }

        // compute the decryption key
        let mut selector: Vec<bool> = Vec::new();
        for _ in 0..t+1 {
            selector.push(true);
        }
        for _ in t+1..n {
            selector.push(false);
        }

        let dec_key = agg_dec(&partial_decryptions, &ct, &selector, &agg_key, &params);
    }
}