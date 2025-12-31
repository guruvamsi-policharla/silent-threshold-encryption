use blstrs::Scalar;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ff::Field;
use silent_threshold_encryption::{
    kzg::KZG10,
    setup::{LagrangePowers, SecretKey},
};

fn bench_setup(c: &mut Criterion) {
    use rand_core::OsRng;
    // WARNING: This benchmark will take a very long time. It is only meant to measure the speedup when compared to the faster Lagrange setup
    let mut group = c.benchmark_group("setup");
    group.sample_size(10);
    let mut rng = OsRng;
    for size in 3..=7 {
        let n = 1 << size; // actually n-1 total parties. one party is a dummy party that is always true
        let tau = Scalar::random(&mut rng);
        let params = KZG10::setup(n, tau).unwrap();

        let sk = SecretKey::new(&mut rng);

        group.bench_with_input(BenchmarkId::from_parameter(n), &params, |b, inp| {
            b.iter(|| sk.get_pk(0, &inp, n));
        });
    }

    group.finish();

    let mut group = c.benchmark_group("Lagrange setup");
    group.sample_size(10);
    let mut rng = OsRng;
    for size in 3..=10 {
        let n = 1 << size; // actually n-1 total parties. one party is a dummy party that is always true
        let tau = Scalar::random(&mut rng);
        let lagrange_params = LagrangePowers::new(tau, n);

        let sk = SecretKey::new(&mut rng);

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
