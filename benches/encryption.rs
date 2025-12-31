use blstrs::Scalar;
use criterion::{criterion_group, criterion_main, Criterion};
use ff::Field;
use silent_threshold_encryption::{
    encryption::encrypt,
    kzg::KZG10,
    setup::{AggregateKey, PublicKey, SecretKey},
};

fn bench_encrypt(c: &mut Criterion) {
    use rand_core::OsRng;
    let mut rng = OsRng;
    let n = 8;
    let t = 2;
    let tau = Scalar::random(&mut rng);
    let params = KZG10::setup(n, tau).unwrap();

    let mut sk: Vec<SecretKey> = Vec::new();
    let mut pk: Vec<PublicKey> = Vec::new();

    for i in 0..n {
        sk.push(SecretKey::new(&mut rng));
        pk.push(sk[i].get_pk(0, &params, n))
    }

    let ak = AggregateKey::new(pk, &params);

    c.bench_function("encrypt", |b| b.iter(|| encrypt(&ak, t, &params)));
}

criterion_group!(benches, bench_encrypt);
criterion_main!(benches);
