use ark_ec::pairing::Pairing;
use ark_poly::univariate::DensePolynomial;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use silent_threshold::{
    encryption::encrypt,
    kzg::KZG10,
    setup::{AggregateKey, PublicKey, SecretKey},
};

type E = ark_bls12_381::Bls12_381;
type UniPoly381 = DensePolynomial<<E as Pairing>::ScalarField>;

//todo: use seeded randomness
fn bench_encrypt(c: &mut Criterion) {
    let mut rng = ark_std::test_rng();
    let n = 8;
    let t = 2;
    let params = KZG10::<E, UniPoly381>::setup(n, &mut rng).unwrap();

    let mut sk: Vec<SecretKey<E>> = Vec::new();
    let mut pk: Vec<PublicKey<E>> = Vec::new();

    for i in 0..n {
        sk.push(SecretKey::<E>::new(&mut rng));
        pk.push(sk[i].get_pk(0, &params, n))
    }

    let ak = AggregateKey::<E>::new(pk, &params);
    // let ct = encrypt::<E>(&ak, t, &params);

    let mut group = c.benchmark_group("encrypt");
    group.bench_with_input(
        BenchmarkId::from_parameter(n),
        &(ak, t, params),
        |b, inp| {
            b.iter(|| encrypt::<E>(&inp.0, inp.1, &inp.2));
        },
    );

    group.finish();
}

criterion_group!(benches, bench_encrypt);
criterion_main!(benches);
