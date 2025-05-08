use ark_ec::pairing::Pairing;
use ark_std::{end_timer, start_timer, Zero};
use silent_threshold_encryption::{
    crs::CRS,
    decryption::agg_dec,
    encryption::encrypt,
    setup::{AggregateKey, SecretKey},
};

type E = ark_bls12_381::Bls12_381;
type G2 = <E as Pairing>::G2;
use rand::seq::IteratorRandom;

fn main() {
    let mut rng = ark_std::test_rng();
    let n = 1 << 4; // actually n-1 total parties. one party is a dummy party that is always true
    let t: usize = 9;
    debug_assert!(t < n);

    let kzg_timer = start_timer!(|| "Setting up KZG parameters");
    let crs = CRS::new(n, &mut rng);
    end_timer!(kzg_timer);

    println!("Setting up key pairs for {} parties", n);
    let key_timer = start_timer!(|| "Setting up keys");

    let sk = (0..n)
        .map(|_| SecretKey::<E>::new(&mut rng))
        .collect::<Vec<_>>();
    let pk = sk
        .iter()
        .enumerate()
        .map(|(i, sk)| sk.get_lagrange_pk(i, &crs))
        .collect::<Vec<_>>();

    end_timer!(key_timer);

    let agg_key_timer = start_timer!(|| "Computing the aggregate key");
    let agg_key = AggregateKey::<E>::new(pk, &crs);
    end_timer!(agg_key_timer);

    let msg = b"Hello, world!";

    let enc_timer = start_timer!(|| "Encrypting a message");
    let ct = encrypt::<E>(&agg_key, t, &crs, msg);
    end_timer!(enc_timer);

    println!("Computing partial decryptions");
    // sample t random signers from 0..n-1
    let mut rng = rand::rng(); // Create a random number generator
    let signers = (0..n).choose_multiple(&mut rng, t);

    let mut selector: Vec<bool> = vec![false; n];
    let mut partial_decryptions: Vec<G2> = vec![G2::zero(); n];

    for i in signers {
        selector[i] = true;
        partial_decryptions[i] = sk[i].partial_decryption(&ct);
    }

    let dec_timer = start_timer!(|| "Aggregating partial decryptions and decrypting");
    let dec_key = agg_dec(&partial_decryptions, &ct, &selector, &agg_key, &crs);
    end_timer!(dec_timer);
    assert_eq!(dec_key, msg, "Decryption failed!");
    println!("Decryption successful!");
}
