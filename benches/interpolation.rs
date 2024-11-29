use ark_ec::pairing::Pairing;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_std::One;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use silent_threshold_encryption::utils::interp_mostly_zero;

type E = ark_bls12_381::Bls12_381;
type F = <E as Pairing>::ScalarField;

fn bench_interpolate(c: &mut Criterion) {
    let mut group = c.benchmark_group("interpolate");

    for size in 3..=10 {
        let n = 1 << size; // actually n-1 total parties. one party is a dummy party that is always true
        let t: usize = n / 2;

        let mut selector: Vec<bool> = Vec::new();
        for _ in 0..t + 1 {
            selector.push(true);
        }
        for _ in t + 1..n {
            selector.push(false);
        }

        let domain = Radix2EvaluationDomain::<F>::new(n).unwrap();
        let domain_elements: Vec<F> = domain.elements().collect();

        let mut points = vec![domain_elements[0]]; // 0 is the dummy party that is always true
        let mut parties: Vec<usize> = Vec::new(); // parties indexed from 0..n-1
        for i in 0..n {
            if selector[i] {
                parties.push(i);
            } else {
                points.push(domain_elements[i]);
            }
        }

        // compute the decryption key
        group.bench_with_input(BenchmarkId::from_parameter(n), &points, |b, inp| {
            b.iter(|| interp_mostly_zero(F::one(), &inp));
        });
    }

    group.finish();
}

criterion_group!(benches, bench_interpolate);
criterion_main!(benches);
