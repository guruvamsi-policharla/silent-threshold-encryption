use ark_std::{end_timer, start_timer};
use blstrs::{G2Projective, Scalar};
use ff::Field;
use group::Group;
use rand::seq::IteratorRandom;
use rayon::prelude::*;
use silent_threshold_encryption::{
    decryption::agg_dec,
    encryption::encrypt,
    kzg::KZG10,
    setup::{AggregateKey, LagrangePowers, SecretKey},
};

fn main() {
    use rand_core::OsRng;
    let mut rng = OsRng;

    let n = 1 << 11; // actually n-1 total parties. one party is a dummy party that is always true
    let t: usize = 1400;
    debug_assert!(t < n);

    let kzg_timer = start_timer!(|| "Setting up KZG parameters");
    let tau = Scalar::random(&mut rng);
    let kzg_params = KZG10::setup(n, tau).unwrap();
    end_timer!(kzg_timer);

    let lagrange_params_timer = start_timer!(|| "Preprocessing lagrange powers");
    let lagrange_params = LagrangePowers::new(tau, n);
    end_timer!(lagrange_params_timer);

    println!("Setting up key pairs for {} parties", n);
    let key_timer = start_timer!(|| "Setting up keys");
    // create the dummy party's keys
    let mut sk = (0..n).map(|_| SecretKey::new(&mut rng)).collect::<Vec<_>>();
    sk[0].nullify();
    let mut pk = vec![sk[0].lagrange_get_pk(0, &lagrange_params, n); n];

    pk.par_iter_mut().enumerate().for_each(|(i, pk_i)| {
        if i > 0 {
            *pk_i = sk[i].lagrange_get_pk(i, &lagrange_params, n);
        }
    });
    end_timer!(key_timer);

    let agg_key_timer = start_timer!(|| "Computing the aggregate key");
    let agg_key = AggregateKey::new(pk, &kzg_params);
    end_timer!(agg_key_timer);

    let enc_timer = start_timer!(|| "Encrypting a message");
    let ct = encrypt(&agg_key, t, &kzg_params);
    end_timer!(enc_timer);

    println!("Computing partial decryptions");
    // sample t random signers from 1..n
    let mut rng = rand::rng(); // Create a random number generator
    let signers = (1..n).choose_multiple(&mut rng, t);

    let mut selector: Vec<bool> = vec![false; n];
    let mut partial_decryptions: Vec<G2Projective> = vec![G2Projective::identity(); n];
    selector[0] = true;
    partial_decryptions[0] = sk[0].partial_decryption(&ct);
    for i in signers {
        selector[i] = true;
        partial_decryptions[i] = sk[i].partial_decryption(&ct);
    }

    let dec_timer = start_timer!(|| "Aggregating partial decryptions and decrypting");
    let dec_key = agg_dec(&partial_decryptions, &ct, &selector, &agg_key, &kzg_params);
    end_timer!(dec_timer);
    assert_eq!(dec_key, ct.enc_key, "Decryption failed!");
    println!("Decryption successful!");
}
