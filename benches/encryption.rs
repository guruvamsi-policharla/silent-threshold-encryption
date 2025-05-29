use criterion::{criterion_group, criterion_main, Criterion};
use silent_threshold_encryption::{
    aggregate::AggregateKey, crs::CRS, encryption::encrypt, setup::SecretKey,
};

type E = ark_bls12_381::Bls12_381;
type G2 = <E as ark_ec::pairing::Pairing>::G2;
use ark_std::UniformRand;

fn bench_encrypt(c: &mut Criterion) {
    let mut rng = ark_std::test_rng();
    let n = 8;
    let t = 2;
    let crs = CRS::<E>::new(n, &mut rng);

    let sk = (0..n)
        .map(|i| SecretKey::<E>::new(&mut rng, i))
        .collect::<Vec<_>>();

    let pk = sk
        .iter()
        .enumerate()
        .map(|(i, sk)| sk.get_lagrange_pk(i, &crs))
        .collect::<Vec<_>>();

    let (_ak, ek) = AggregateKey::<E>::new(pk, &crs);
    let msg = b"Hello, world!";

    let gamma_g2 = G2::rand(&mut rng);

    c.bench_function("encrypt", |b| {
        b.iter(|| encrypt::<E>(&ek, t, &crs, gamma_g2, msg))
    });
}

criterion_group!(benches, bench_encrypt);
criterion_main!(benches);
