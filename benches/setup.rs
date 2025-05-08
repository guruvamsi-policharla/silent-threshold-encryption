use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use silent_threshold_encryption::{crs::CRS, setup::SecretKey};

type E = ark_bls12_381::Bls12_381;

fn bench_setup(c: &mut Criterion) {
    // WARNING: This benchmark will take a very long time. It is only meant to measure the speedup when compared to the faster Lagrange setup
    let mut group = c.benchmark_group("setup");
    group.sample_size(10);
    let mut rng = ark_std::test_rng();
    for size in 3..=7 {
        let n = 1 << size;
        let crs = CRS::<E>::new(n, &mut rng);
        let sk = SecretKey::<E>::new(&mut rng);

        group.bench_with_input(BenchmarkId::from_parameter(n), &crs, |b, inp| {
            b.iter(|| sk.get_lagrange_pk(0, &inp));
        });
    }

    group.finish();
}

criterion_group!(benches, bench_setup);
criterion_main!(benches);
