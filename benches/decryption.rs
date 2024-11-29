use ark_ec::pairing::Pairing;
use ark_poly::univariate::DensePolynomial;
use ark_std::{UniformRand, Zero};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use silent_threshold_encryption::{
    decryption::agg_dec,
    encryption::encrypt,
    kzg::KZG10,
    setup::{AggregateKey, LagrangePowers, PublicKey, SecretKey},
};

type E = ark_bls12_381::Bls12_381;
type G2 = <E as Pairing>::G2;
type Fr = <E as Pairing>::ScalarField;
type UniPoly381 = DensePolynomial<<E as Pairing>::ScalarField>;

fn bench_decrypt(c: &mut Criterion) {
    let mut rng = ark_std::test_rng();
    let mut group = c.benchmark_group("decrypt");

    for size in 3..=10 {
        let n = 1 << size; // actually n-1 total parties. one party is a dummy party that is always true
        let t: usize = n / 2;

        let tau = Fr::rand(&mut rng);
        let params = KZG10::<E, UniPoly381>::setup(n, tau.clone()).unwrap();
        let lagrange_params = LagrangePowers::<E>::new(tau, n);

        let mut sk: Vec<SecretKey<E>> = Vec::new();
        let mut pk: Vec<PublicKey<E>> = Vec::new();

        // create the dummy party's keys
        sk.push(SecretKey::<E>::new(&mut rng));
        sk[0].nullify();
        pk.push(sk[0].lagrange_get_pk(0, &lagrange_params, n));

        for i in 1..n {
            sk.push(SecretKey::<E>::new(&mut rng));
            pk.push(sk[i].lagrange_get_pk(i, &lagrange_params, n));
        }

        let agg_key = AggregateKey::<E>::new(pk, &params);
        let ct = encrypt::<E>(&agg_key, t, &params);

        // compute partial decryptions
        let mut partial_decryptions: Vec<G2> = Vec::new();
        for i in 0..t + 1 {
            partial_decryptions.push(sk[i].partial_decryption(&ct));
        }
        for _ in t + 1..n {
            partial_decryptions.push(G2::zero());
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
