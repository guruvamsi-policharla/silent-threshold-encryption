use criterion::{criterion_group, criterion_main, Criterion};
use silent_threshold_encryption::{
    crs::CRS,
    encryption::encrypt,
    setup::{AggregateKey, SecretKey},
};

type E = ark_bls12_381::Bls12_381;

fn bench_encrypt(c: &mut Criterion) {
    let mut rng = ark_std::test_rng();
    let n = 8;
    let t = 2;
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

    c.bench_function("encrypt", |b| {
        b.iter(|| encrypt::<E>(&agg_key, t, &crs, msg))
    });
}

criterion_group!(benches, bench_encrypt);
criterion_main!(benches);
