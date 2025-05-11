use ark_ec::pairing::Pairing;
use ark_std::{end_timer, start_timer, Zero};
use silent_threshold_encryption::{
    aggregate::SystemPublicKeys,
    crs::CRS,
    decryption::agg_dec,
    encryption::encrypt,
    setup::{LagPolys, SecretKey},
};

type E = ark_bls12_381::Bls12_381;
type G2 = <E as Pairing>::G2;
use rand::seq::IteratorRandom;

fn main() {
    let mut rng = ark_std::test_rng();
    let n = 1 << 3;
    let m = n * (1 << 0);
    let t: usize = n / 2;
    debug_assert!(t < n);
    let k = 3;

    let kzg_timer = start_timer!(|| "Setting up parameters");
    let crs = CRS::new(n, &mut rng);
    let lag_polys = LagPolys::new(n);
    end_timer!(kzg_timer);

    println!("Setting up key pairs for {} parties", m);
    let sk = (0..m)
        .map(|_| SecretKey::<E>::new(&mut rng))
        .collect::<Vec<_>>();
    let pk = sk
        .iter()
        .enumerate()
        .map(|(i, sk)| sk.get_pk(i, &crs))
        .collect::<Vec<_>>();

    let setup_timer = start_timer!(|| "Setting up system keys");
    let system_keys = SystemPublicKeys::<E>::new(pk.clone(), &crs, &lag_polys, k);
    end_timer!(setup_timer);

    let subset_timer = start_timer!(|| "Computing the aggregate key of a subset");
    let subset = (0..n).collect::<Vec<_>>();
    let subset_agg_key = system_keys.get_aggregate_key(&subset, &crs, &lag_polys);
    end_timer!(subset_timer);

    let msg = b"Hello, world!";

    let enc_timer = start_timer!(|| "Encrypting a message");
    let ct = encrypt::<E>(&subset_agg_key, t, &crs, msg);
    end_timer!(enc_timer);

    println!("Computing partial decryptions");
    // // sample t random signers from 0..n-1
    // let mut rng = rand::rng(); // Create a random number generator
    // let signers = (0..n).choose_multiple(&mut rng, t);

    // let mut selector: Vec<bool> = vec![false; n];
    // let mut partial_decryptions: Vec<G2> = vec![G2::zero(); n];

    // for i in signers {
    //     selector[i] = true;
    //     partial_decryptions[i] = sk[i].partial_decryption(&ct);
    // }

    // compute partial decryptions
    let mut partial_decryptions: Vec<G2> = Vec::new();
    for i in 0..t {
        let id = subset_agg_key.lag_pks[i].id;
        partial_decryptions.push(sk[id].partial_decryption(&ct));
    }
    for _ in t..n {
        partial_decryptions.push(G2::zero());
    }

    // compute the decryption key
    let mut selector: Vec<bool> = Vec::new();
    for _ in 0..t {
        selector.push(true);
    }
    for _ in t..n {
        selector.push(false);
    }

    let dec_timer = start_timer!(|| "Aggregating partial decryptions and decrypting");
    let dec_key = agg_dec(&partial_decryptions, &ct, &selector, &subset_agg_key, &crs);
    end_timer!(dec_timer);
    assert_eq!(dec_key, msg, "Decryption failed!");
    println!("Decryption successful!");
}
