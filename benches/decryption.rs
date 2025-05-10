use ark_ec::pairing::Pairing;
use ark_std::Zero;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use silent_threshold_encryption::{
    aggregate::AggregateKey, crs::CRS, decryption::agg_dec, encryption::encrypt, setup::SecretKey,
};

type E = ark_bls12_381::Bls12_381;
type G2 = <E as Pairing>::G2;

fn bench_decrypt(c: &mut Criterion) {
    let mut rng = ark_std::test_rng();
    let mut group = c.benchmark_group("decrypt");

    for size in 3..=10 {
        let n = 1 << size; // actually n-1 total parties. one party is a dummy party that is always true
        let t: usize = n / 2;

        let crs = CRS::<E>::new(n, &mut rng);

        let sk = (0..n)
            .map(|_| SecretKey::<E>::new(&mut rng))
            .collect::<Vec<_>>();

        let pk = sk
            .iter()
            .enumerate()
            .map(|(i, sk)| sk.get_lagrange_pk(i, &crs))
            .collect::<Vec<_>>();

        let agg_key = AggregateKey::<E>::new(pk, &crs);
        let msg = b"Hello, world!";
        let ct = encrypt::<E>(&agg_key, t, &crs, msg);

        // compute partial decryptions
        let mut partial_decryptions: Vec<G2> = Vec::new();
        for i in 0..t {
            partial_decryptions.push(sk[i].partial_decryption(&ct));
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

        group.bench_with_input(
            BenchmarkId::from_parameter(n),
            &(partial_decryptions, ct, selector, agg_key, crs),
            |b, inp| {
                b.iter(|| agg_dec(&inp.0, &inp.1, &inp.2, &inp.3, &inp.4));
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_decrypt);
criterion_main!(benches);
