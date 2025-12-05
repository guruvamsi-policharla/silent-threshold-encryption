use blstrs::{G1Affine, G1Projective, G2Affine, G2Prepared, G2Projective, Gt, Scalar};
use ff::Field;
use group::{Curve, Group};
use pairing::{MillerLoopResult, MultiMillerLoop};

use crate::{
    encryption::Ciphertext,
    kzg::{PowersOfTau, KZG10},
    polynomial::{DensePolynomial, Radix2EvaluationDomain},
    setup::AggregateKey,
    utils::interp_mostly_zero,
};

pub fn agg_dec(
    partial_decryptions: &[G2Projective], //insert 0 if a party did not respond or verification failed
    ct: &Ciphertext,
    selector: &[bool],
    agg_key: &AggregateKey,
    params: &PowersOfTau,
) -> Gt {
    let n = agg_key.pk.len();
    let domain = Radix2EvaluationDomain::new(n).unwrap();
    let domain_elements: Vec<Scalar> = domain.elements().collect();

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

    let b = interp_mostly_zero(Scalar::ONE, &points);
    let b_evals = domain.fft(&b.coeffs);

    debug_assert!(b.degree() == points.len() - 1);
    debug_assert!(b.evaluate(&domain_elements[0]) == Scalar::ONE);

    // commit to b in g2
    let b_g2: G2Projective = KZG10::commit_g2(params, &b).unwrap().into();

    // q0 = (b-1)/(x-domain_elements[0])
    let mut bminus1 = b.clone();
    bminus1.coeffs[0] -= Scalar::ONE;

    debug_assert!(bminus1.evaluate(&domain_elements[0]) == Scalar::ZERO);

    let xminus1 = DensePolynomial::from_coefficients_vec(vec![-domain_elements[0], Scalar::ONE]);
    let q0 = &bminus1 / &xminus1;

    let q0_g1: G1Projective = KZG10::commit_g1(params, &q0).unwrap().into();

    // bhat = x^{t+1} * b
    // insert t+1 0s at the beginning of bhat.coeffs
    let mut bhat_coeffs = vec![Scalar::ZERO; ct.t + 1];
    bhat_coeffs.append(&mut b.coeffs.clone());
    let bhat = DensePolynomial::from_coefficients_vec(bhat_coeffs);
    debug_assert_eq!(bhat.degree(), n);

    let bhat_g1: G1Projective = KZG10::commit_g1(params, &bhat).unwrap().into();

    let n_inv = Scalar::from(n as u64).invert().unwrap();

    // compute the aggregate public key using MSM
    let mut apk = G1Projective::identity();
    for &i in &parties {
        apk += agg_key.pk[i].bls_pk * b_evals[i];
    }
    apk *= n_inv;

    // compute sigma = (\sum B(omega^i)partial_decryptions[i])/(n) for i in parties
    let mut sigma = G2Projective::identity();
    for &i in &parties {
        sigma += partial_decryptions[i] * b_evals[i];
    }
    sigma *= n_inv;

    // compute Qx, Qhatx and Qz
    let mut qx = G1Projective::identity();
    for &i in &parties {
        qx += agg_key.pk[i].sk_li_x * b_evals[i];
    }

    let mut qz = G1Projective::identity();
    for &i in &parties {
        qz += agg_key.agg_sk_li_lj_z[i] * b_evals[i];
    }

    let mut qhatx = G1Projective::identity();
    for &i in &parties {
        qhatx += agg_key.pk[i].sk_li_minus0 * b_evals[i];
    }

    // e(w1||sa1, sa2||w2)
    let minus1 = -Scalar::ONE;
    let w1 = [
        apk * minus1,
        qz * minus1,
        qx * minus1,
        qhatx,
        bhat_g1 * minus1,
        q0_g1 * minus1,
    ];
    let w2 = [b_g2, sigma];

    let mut enc_key_lhs = w1.to_vec();
    enc_key_lhs.extend_from_slice(&ct.sa1);

    let mut enc_key_rhs = ct.sa2.to_vec();
    enc_key_rhs.extend_from_slice(&w2);

    // Convert to affine for pairing
    let lhs_affine: Vec<G1Affine> = enc_key_lhs.iter().map(|p| p.to_affine()).collect();
    let rhs_affine: Vec<G2Affine> = enc_key_rhs.iter().map(|p| p.to_affine()).collect();

    // Prepare G2 elements for pairing
    let rhs_prepared: Vec<G2Prepared> = rhs_affine.iter().map(|p| G2Prepared::from(*p)).collect();

    // Compute multi-pairing
    let terms: Vec<(&G1Affine, &G2Prepared)> = lhs_affine.iter().zip(rhs_prepared.iter()).collect();

    let enc_key = blstrs::Bls12::multi_miller_loop(&terms).final_exponentiation();

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
    use rand_core::OsRng;

    #[test]
    fn test_decryption() {
        let mut rng = OsRng;
        let n = 1 << 4; // actually n-1 total parties. one party is a dummy party that is always true
        let t: usize = n / 2;
        debug_assert!(t < n);

        let tau = Scalar::random(&mut rng);
        let params = KZG10::setup(n, tau).unwrap();

        let mut sk: Vec<SecretKey> = Vec::new();
        let mut pk: Vec<PublicKey> = Vec::new();

        // create the dummy party's keys
        sk.push(SecretKey::new(&mut rng));
        sk[0].nullify();
        pk.push(sk[0].get_pk(0, &params, n));

        for i in 1..n {
            sk.push(SecretKey::new(&mut rng));
            pk.push(sk[i].get_pk(i, &params, n))
        }

        let agg_key = AggregateKey::new(pk, &params);
        let ct = encrypt(&agg_key, t, &params);

        // compute partial decryptions
        let mut partial_decryptions: Vec<G2Projective> = Vec::new();
        for i in 0..t + 1 {
            partial_decryptions.push(sk[i].partial_decryption(&ct));
        }
        for _ in t + 1..n {
            partial_decryptions.push(G2Projective::identity());
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
