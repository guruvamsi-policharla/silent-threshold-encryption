use ark_ec::pairing::Pairing;
use ark_poly::univariate::DensePolynomial;
use ark_std::{UniformRand, Zero};
use silent_threshold::{
    decryption::agg_dec,
    encryption::encrypt,
    kzg::KZG10,
    setup::{AggregateKey, LagrangePowers, PublicKey, SecretKey},
};

type E = ark_bls12_381::Bls12_381;
type G2 = <E as Pairing>::G2;
type Fr = <E as Pairing>::ScalarField;
type UniPoly381 = DensePolynomial<<E as Pairing>::ScalarField>;

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
    // compute partial decryptions
    let mut partial_decryptions: Vec<G2> = Vec::new();
    for i in 0..t + 1 {
        partial_decryptions.push(sk[i].partial_decryption(&ct));
    }
    for _ in t + 1..n {
        partial_decryptions.push(G2::zero());
    }

    println!("Aggregating partial decryptions and decrypting");
    // compute the decryption key
    let mut selector: Vec<bool> = Vec::new();
    for _ in 0..t + 1 {
        selector.push(true);
    }
    for _ in t + 1..n {
        selector.push(false);
    }

    let _dec_key = agg_dec(&partial_decryptions, &ct, &selector, &agg_key, &kzg_params);

    println!("Decryption successful!");
}
