use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use silent_threshold_encryption::{
    aggregate::AggregateKey,
    crs::CRS,
    decryption::agg_dec,
    encryption::encrypt,
    setup::{PartialDecryption, SecretKey},
};

type E = ark_bls12_381::Bls12_381;
type G2 = <E as ark_ec::pairing::Pairing>::G2;
use ark_std::UniformRand;

fn bench_decrypt(c: &mut Criterion) {
    let mut rng = ark_std::test_rng();
    let mut group = c.benchmark_group("decrypt");

    for size in 3..=10 {
        let n = 1 << size; // actually n-1 total parties. one party is a dummy party that is always true
        let t: usize = n / 2;

        let crs = CRS::<E>::new(n, &mut rng);

        let sk = (0..n)
            .map(|i| SecretKey::<E>::new(&mut rng, i))
            .collect::<Vec<_>>();

        let pk = sk
            .iter()
            .enumerate()
            .map(|(i, sk)| sk.get_lagrange_pk(i, &crs))
            .collect::<Vec<_>>();

        let (ak, ek) = AggregateKey::<E>::new(pk, &crs);
        let msg = b"Hello, world!";
        let gamma_g2 = G2::rand(&mut rng);
        let ct = encrypt::<E>(&ek, t, &crs, gamma_g2, msg);

        // compute partial decryptions
        let mut partial_decryptions: Vec<PartialDecryption<E>> = Vec::new();
        for i in 0..t {
            partial_decryptions.push(sk[i].partial_decryption(&ct));
        }
        for _ in t..n {
            partial_decryptions.push(PartialDecryption::zero());
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
            &(partial_decryptions, ct, selector, ak, crs),
            |b, inp| {
                b.iter(|| agg_dec(&inp.0, &inp.1, &inp.2, &inp.3, &inp.4));
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_decrypt);
criterion_main!(benches);
