use blstrs::{G2Projective, Scalar};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ff::Field;
use group::Group;
use silent_threshold_encryption::{
    decryption::agg_dec,
    encryption::encrypt,
    kzg::KZG10,
    setup::{AggregateKey, LagrangePowers, PublicKey, SecretKey},
};

fn bench_decrypt(c: &mut Criterion) {
    use rand_core::OsRng;
    let mut rng = OsRng;
    let mut group = c.benchmark_group("decrypt");

    for size in 3..=10 {
        let n = 1 << size; // actually n-1 total parties. one party is a dummy party that is always true
        let t: usize = n / 2;

        let tau = Scalar::random(&mut rng);
        let params = KZG10::setup(n, tau).unwrap();
        let lagrange_params = LagrangePowers::new(tau, n);

        let mut sk: Vec<SecretKey> = Vec::new();
        let mut pk: Vec<PublicKey> = Vec::new();

        // create the dummy party's keys
        sk.push(SecretKey::new(&mut rng));
        sk[0].nullify();
        pk.push(sk[0].lagrange_get_pk(0, &lagrange_params, n));

        for i in 1..n {
            sk.push(SecretKey::new(&mut rng));
            pk.push(sk[i].lagrange_get_pk(i, &lagrange_params, n));
        }

        let agg_key = AggregateKey::new(pk, &params);
        let ct = encrypt(&agg_key, t, &params);

        // compute partial decryptions
        let mut partial_decryptions: Vec<G2Projective> = Vec::new();
        for i in 0..t + 1 {
            partial_decryptions.push(sk[i].partial_decryption(&ct));
        }
        for _ in t + 1..n {
            partial_decryptions.push(G2Projective::identity());
        }

        // compute the decryption key
        let mut selector: Vec<bool> = Vec::new();
        for _ in 0..t + 1 {
            selector.push(true);
        }
        for _ in t + 1..n {
            selector.push(false);
        }

        group.bench_with_input(
            BenchmarkId::from_parameter(n),
            &(partial_decryptions, ct, selector, agg_key, params),
            |b, inp| {
                b.iter(|| agg_dec(&inp.0, &inp.1, &inp.2, &inp.3, &inp.4));
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_decrypt);
criterion_main!(benches);
