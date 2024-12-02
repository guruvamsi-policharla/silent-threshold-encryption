use ark_ec::pairing::Pairing;
use ark_poly::univariate::DensePolynomial;
use ark_std::{UniformRand, Zero};
use silent_threshold_encryption::{
    decryption::agg_dec,
    encryption::encrypt,
    kzg::KZG10,
    setup::{AggregateKey, LagrangePowers, PublicKey, SecretKey},
};

type E = ark_bls12_381::Bls12_381;
type G2 = <E as Pairing>::G2;
type Fr = <E as Pairing>::ScalarField;
type UniPoly381 = DensePolynomial<<E as Pairing>::ScalarField>;
use rand::{seq::IteratorRandom, thread_rng};

fn main() {
    let mut rng = ark_std::test_rng();
    let n = 1 << 5; // actually n-1 total parties. one party is a dummy party that is always true
    let t: usize = 9;
    debug_assert!(t < n);

    println!("Setting up KZG parameters");
    let tau = Fr::rand(&mut rng);
    let kzg_params = KZG10::<E, UniPoly381>::setup(n, tau.clone()).unwrap();

    println!("Preprocessing lagrange powers");
    let lagrange_params = LagrangePowers::<E>::new(tau, n);

    println!("Setting up key pairs for {} parties", n);
    let mut sk: Vec<SecretKey<E>> = Vec::new();
    let mut pk: Vec<PublicKey<E>> = Vec::new();

    // create the dummy party's keys
    sk.push(SecretKey::<E>::new(&mut rng));
    sk[0].nullify();
    pk.push(sk[0].lagrange_get_pk(0, &lagrange_params, n));

    for i in 1..n {
        sk.push(SecretKey::<E>::new(&mut rng));
        pk.push(sk[i].lagrange_get_pk(i, &lagrange_params, n))
    }

    println!("Compting the aggregate key");
    let agg_key = AggregateKey::<E>::new(pk, &kzg_params);

    println!("Encrypting a message");
    let ct = encrypt::<E>(&agg_key, t, &kzg_params);

    println!("Computing partial decryptions");

    // sample t random signers from 1..n
    let mut rng = thread_rng(); // Create a random number generator
    let signers = (1..n).choose_multiple(&mut rng, t);

    let mut selector: Vec<bool> = vec![false; n];
    let mut partial_decryptions: Vec<G2> = vec![G2::zero(); n];
    selector[0] = true;
    partial_decryptions[0] = sk[0].partial_decryption(&ct);
    for i in signers {
        selector[i] = true;
        partial_decryptions[i] = sk[i].partial_decryption(&ct);
    }

    println!("Aggregating partial decryptions and decrypting");
    let dec_key = agg_dec(&partial_decryptions, &ct, &selector, &agg_key, &kzg_params);
    assert_eq!(dec_key, ct.enc_key, "Decryption failed!");
    println!("Decryption successful!");
}
