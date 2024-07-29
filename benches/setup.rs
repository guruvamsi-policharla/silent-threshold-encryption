use ark_ec::pairing::Pairing;
use ark_poly::univariate::DensePolynomial;
use ark_std::UniformRand;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use silent_threshold::{
    kzg::KZG10,
    setup::{LagrangePowers, SecretKey},
};

type E = ark_bls12_381::Bls12_381;
type Fr = <E as Pairing>::ScalarField;
type UniPoly381 = DensePolynomial<<E as Pairing>::ScalarField>;

fn bench_setup(c: &mut Criterion) {
    // WARNING: This benchmark will take a very long time. It is only meant to measure the speedup when compared to the faster Lagrange setup
    let mut group = c.benchmark_group("setup");
    group.sample_size(10);
    let mut rng = ark_std::test_rng();
    for size in 3..=7 {
        let n = 1 << size; // actually n-1 total parties. one party is a dummy party that is always true
        let tau = Fr::rand(&mut rng);
        let params = KZG10::<E, UniPoly381>::setup(n, tau.clone()).unwrap();

        let sk = SecretKey::<E>::new(&mut rng);

        group.bench_with_input(BenchmarkId::from_parameter(n), &params, |b, inp| {
            b.iter(|| sk.get_pk(0, &inp, n));
        });
    }

    group.finish();

    let mut group = c.benchmark_group("Lagrange setup");
    group.sample_size(10);
    let mut rng = ark_std::test_rng();
    for size in 3..=10 {
        let n = 1 << size; // actually n-1 total parties. one party is a dummy party that is always true
        let tau = Fr::rand(&mut rng);
        let lagrange_params = LagrangePowers::<E>::new(tau, n);

        let sk = SecretKey::<E>::new(&mut rng);

        group.bench_with_input(
            BenchmarkId::from_parameter(n),
            &lagrange_params,
            |b, inp| {
                b.iter(|| sk.lagrange_get_pk(0, &inp, n));
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_setup);
criterion_main!(benches);
