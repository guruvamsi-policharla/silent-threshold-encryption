use std::ops::{Sub, Div, Mul};

use ark_ec::pairing::{Pairing, PairingOutput};
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial,
    Radix2EvaluationDomain,
};
use ark_std::{One, Zero};

use crate::{
    encryption::Ciphertext,
    kzg::{UniversalParams, KZG10},
    setup::AggregateKey,
    utils::{interp_mostly_zero, lagrange_poly},
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
    let mut points = vec![domain_elements[0]]; // 0 is the dummy party that is always true
    for i in 0..n {
        if selector[i] == false {
            points.push(domain_elements[i]);
            println!("i: {}", i);
        }
    }

    let b = interp_mostly_zero(E::ScalarField::one(), &points);
    let b_evals = domain.fft(&b.coeffs);

    debug_assert!(b.degree() == points.len() - 1);
    // println!("b.degree(): {}", b.degree());
    // println!("points: {:?}", points);
    debug_assert!(b.evaluate(&domain_elements[0]) == E::ScalarField::one());

    // commit to b in g2
    let b_g2: E::G2 = KZG10::<E, DensePolynomial<E::ScalarField>>::commit_g2(&params, &b)
        .unwrap()
        .into();

    // q0 = (b-1)/(x-domain_elements[0])
    let mut bminus1 = b.clone();
    bminus1.coeffs[0] -= E::ScalarField::one();

    debug_assert!(bminus1.evaluate(&domain_elements[0]) == E::ScalarField::zero());

    let xminus1 =
        DensePolynomial::from_coefficients_vec(vec![-domain_elements[0], E::ScalarField::one()]);
    let q0 = bminus1.div(&xminus1);

    let q0_g1: E::G1 = KZG10::<E, DensePolynomial<E::ScalarField>>::commit_g1(&params, &q0)
        .unwrap()
        .into();

    // check that KZG verification passes
    // let htau_x = params.powers_of_h[1] + (params.powers_of_h[0] * (-domain.element(0)));
    // let rhs = KZG10::commit_g2(&params, &xminus1).unwrap().into();
    // debug_assert_eq!(htau_x, rhs);

    // debug_assert_eq!(E::pairing(E::G1::generator(), E::G2::generator()), E::pairing(params.powers_of_g[0], params.powers_of_h[0]));

    // debug_assert_eq!(
    //     E::pairing(
    //         params.powers_of_g[0],
    //         b_g2 + (params.powers_of_h[0]*(-E::ScalarField::one())),
    //     ),
    //     E::pairing(q0_g1, htau_x)
    // );

    // bhat = x^t * b
    // insert t 0s at the beginning of bhat.coeffs
    let mut bhat_coeffs = vec![E::ScalarField::zero(); ct.t];
    bhat_coeffs.append(&mut b.coeffs.clone());
    let bhat = DensePolynomial::from_coefficients_vec(bhat_coeffs);
    debug_assert_eq!(bhat.degree(), n - 1);
    // println!("bhat.degree():{}", bhat.degree());

    let bhat_g1: E::G1 = KZG10::<E, DensePolynomial<E::ScalarField>>::commit_g1(&params, &bhat)
        .unwrap()
        .into();

    let n_inv = E::ScalarField::one() / E::ScalarField::from((n) as u32);
    // compute the aggregate public key
    let mut apk: E::G1 = E::G1::zero();
    for i in 0..agg_key.pk.len() {
        if selector[i] {
            apk += agg_key.pk[i].bls_pk * b_evals[i];
        }
    }
    apk *= n_inv;

    // compute sigma = (\sum B(omega^i)partial_decryptions[i])/(n) for i in parties
    let mut sigma: E::G2 = E::G2::zero();
    for i in 0..agg_key.pk.len() {
        if selector[i] {
            sigma += partial_decryptions[i] * b_evals[i];
        }
    }
    sigma *= n_inv;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // todo: optimize this by populating parties earlier along with points
    // todo: avoid a linear size traversal and only do RAM access
    // compute Qx, Qhatx and Qz
    let mut parties: Vec<usize> = Vec::new(); // parties indexed from 0..n-1
    for i in 0..n {
        if selector[i] {
            parties.push(i);
        }
    }

    println!("parties:{:?}", parties);
    
    let mut qx: E::G1 = E::G1::zero();
    for &i in &parties {
        qx += agg_key.pk[i].sk_li_by_tau * b_evals[i];
    }

    #[cfg(debug_assertions)] {
        if ct.t==0 {
            println!("checking qx");
            let l0: DensePolynomial<E::ScalarField> = lagrange_poly(n, 0);
            let l0_by_tau = DensePolynomial::<E::ScalarField>::from_coefficients_vec(l0.coeffs[1..].to_vec());
            let should_be_qx: E::G1 = KZG10::<E, DensePolynomial<E::ScalarField>>::commit_g1(&params, &l0_by_tau).unwrap().into();
            debug_assert_eq!(should_be_qx, qx);
        }
    }

    let mut qz: E::G1 = E::G1::zero();
    for &i in &parties {
        let mut qz_i = E::G1::zero();
        for j in 0..n {
            qz_i += agg_key.pk[j].sk_li_by_z[i];
        }
        qz += qz_i * b_evals[i];
    }

    #[cfg(debug_assertions)] {
        if ct.t==0 {
            println!("checking qz");
            let l0: DensePolynomial<E::ScalarField> = lagrange_poly(n, 0);
            let mut qz_poly = l0.mul(&l0).sub(&l0);

            for i in 1..n {
                qz_poly = qz_poly + &l0.mul(&lagrange_poly(n, i))*agg_key.pk[i].sk;
            }

            qz_poly = qz_poly.divide_by_vanishing_poly(domain).unwrap().0;

            let should_be_qz: E::G1 = KZG10::<E, DensePolynomial<E::ScalarField>>::commit_g1(&params, &qz_poly).unwrap().into();
            debug_assert_eq!(qz, should_be_qz);
        }
    }

    let mut qhatx: E::G1 = E::G1::zero();
    for &i in &parties {
        qhatx += agg_key.pk[i].sk_li_minus0 * b_evals[i];
    }

    println!("qx:{:?}", qx);
    println!("qz:{:?}", qz);
    println!("qhatx:{:?}", qhatx);

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
    let w2 = [
        b_g2, 
        sigma
    ];

    // e(-q0_g1, sa2[5]) + e(sa1[0], b_g2) = enc_key when s0=s1=s2=s3=0
    // e(-w1[5], sa2[5]) + e(sa1[0], w2[0]) = enc_key when s0=s1=s2=s3=0
    // let lhs = E::pairing(q0_g1*(minus1), ct.sa2[5]);
    // let rhs = E::pairing(ct.sa1[0], b_g2);
    // debug_assert_eq!(lhs+rhs, ct.enc_key);

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

    type E = ark_bls12_381::Bls12_381;
    // type F = <E as Pairing>::ScalarField;
    // type G1 = <E as Pairing>::G1;
    type G2 = <E as Pairing>::G2;
    type UniPoly381 = DensePolynomial<<E as Pairing>::ScalarField>;

    #[test]
    fn test_decryption() {
        let mut rng = ark_std::test_rng();
        let n = 1 << 3; // actually n-1 total parties. one party is a dummy party that is always true
        let t: usize = 6;
        debug_assert!(t < n);

        let params = KZG10::<E, UniPoly381>::setup(n, &mut rng).unwrap();

        let mut sk: Vec<SecretKey<E>> = Vec::new();
        let mut pk: Vec<PublicKey<E>> = Vec::new();

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

        let dec_key = agg_dec(&partial_decryptions, &ct, &selector, &agg_key, &params);
    }
}
