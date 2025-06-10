use ark_std::{end_timer, start_timer};
use silent_threshold_encryption::{
	aggregate::SystemPublicKeys,
	crs::CRS,
	decryption::agg_dec,
	encryption::encrypt,
	setup::{LagPolys, PartialDecryption, SecretKey},
};
type E = ark_bls12_381::Bls12_381;
type G2 = <E as ark_ec::pairing::Pairing>::G2;

use ark_std::UniformRand;
use rand::seq::IteratorRandom;

fn main() {
	let mut rng = ark_std::test_rng();
	let n = 1 << 6;
	let m = n * (1 << 1);
	let t: usize = n / 2;
	debug_assert!(t < n);
	let k = 3;

	let kzg_timer = start_timer!(|| "Setting up parameters");
	let crs = CRS::new(n, &mut rng);
	let lag_polys = LagPolys::new(n);
	end_timer!(kzg_timer);

	println!("Setting up key pairs for {} parties", m);
	let sk = (0..m).map(|i| SecretKey::<E>::new(&mut rng, i)).collect::<Vec<_>>();
	let pk = sk.iter().map(|sk| sk.get_pk(&crs)).collect::<Vec<_>>();

	let setup_timer = start_timer!(|| "Setting up system keys");
	let system_keys = SystemPublicKeys::<E>::new(pk.clone(), &crs, &lag_polys, k);
	end_timer!(setup_timer);

	let subset_timer = start_timer!(|| "Computing the aggregate key of a subset");
	let mut thread_rng = rand::rng(); // Create a random number generator
	let subset = (0..m).choose_multiple(&mut thread_rng, n);
	let (ak, ek) = system_keys.get_aggregate_key(&subset, &crs, &lag_polys);
	end_timer!(subset_timer);

	let msg = b"Hello, world!";

	let enc_timer = start_timer!(|| "Encrypting a message");
	let gamma_g2 = G2::rand(&mut rng);
	let ct = encrypt::<E>(&ek, t, &crs, gamma_g2, msg);
	end_timer!(enc_timer);

	println!("Computing partial decryptions");
	// sample t random signers positions
	let signer_positions = (0..crs.n).choose_multiple(&mut thread_rng, t);

	let mut selector: Vec<bool> = vec![false; n];
	let mut partial_decryptions: Vec<PartialDecryption<E>> = vec![PartialDecryption::zero(); n];

	for i in signer_positions {
		selector[i] = true;
		let id = ak.lag_pks[i].id;
		partial_decryptions[i] = sk[id].partial_decryption(&ct);
	}

	// // compute partial decryptions
	// let mut partial_decryptions: Vec<G2> = Vec::new();
	// for i in 0..t {
	//     let id = ak.lag_pks[i].id;
	//     partial_decryptions.push(sk[id].partial_decryption(&ct));
	// }
	// for _ in t..n {
	//     partial_decryptions.push(G2::zero());
	// }

	// // compute the decryption key
	// let mut selector: Vec<bool> = Vec::new();
	// for _ in 0..t {
	//     selector.push(true);
	// }
	// for _ in t..n {
	//     selector.push(false);
	// }

	let dec_timer = start_timer!(|| "Aggregating partial decryptions and decrypting");
	let dec_key = agg_dec(&partial_decryptions, &ct, &selector, &ak, &crs);
	end_timer!(dec_timer);
	assert_eq!(dec_key, msg, "Decryption failed!");
	println!("Decryption successful!");
}
