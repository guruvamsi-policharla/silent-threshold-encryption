use crate::crs::CRS;
use crate::setup::{LagPolys, LagPublicKey, PublicKey};
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{end_timer, start_timer, One, Zero};
use rand::Rng;
use rand::SeedableRng;

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct AggregateKey<E: Pairing> {
    pub pk: Vec<LagPublicKey<E>>,
    pub agg_sk_li_lj_z: Vec<E::G1>,
    pub ask: E::G1,
    pub z_g2: E::G2,

    //preprocessed values
    pub e_gh: PairingOutput<E>,
}

impl<E: Pairing> AggregateKey<E> {
    pub fn new(pk: Vec<LagPublicKey<E>>, crs: &CRS<E>) -> Self {
        let n = pk.len();
        let z_g2 = crs.powers_of_h[n] + crs.powers_of_h[0] * (-E::ScalarField::one());

        // gather sk_li from all public keys
        let mut ask = E::G1::zero();
        for pki in pk.iter() {
            ask += pki.sk_li;
        }

        let mut agg_sk_li_lj_z = vec![];
        for i in 0..n {
            let mut agg_sk_li_lj_zi = E::G1::zero();
            for pkj in pk.iter() {
                agg_sk_li_lj_zi += pkj.sk_li_lj_z[i];
            }
            agg_sk_li_lj_z.push(agg_sk_li_lj_zi);
        }

        AggregateKey {
            pk,
            agg_sk_li_lj_z,
            ask,
            z_g2,
            e_gh: E::pairing(crs.powers_of_g[0], crs.powers_of_h[0]),
        }
    }
}

/// contains public keys of all m parties in the system
/// also maintains lagrange public keys for each party at k "random" positions
/// this ensures that the public keys of any subset of n parties forms an "almost" perfect matching
/// hence, allowing for efficient aggregation and decryption
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct SystemPublicKeys<E: Pairing> {
    pub m: usize,
    pub k: usize,
    pub pk: Vec<PublicKey<E>>,
    pub lag_pk: Vec<Vec<LagPublicKey<E>>>,
}

impl<E: Pairing> SystemPublicKeys<E> {
    pub fn new(
        pk: Vec<PublicKey<E>>,
        crs: &CRS<E>,
        lag_polys: LagPolys<E::ScalarField>,
        k: usize,
    ) -> Self {
        // using a deterministic seed for reproducibility across machines
        // can derandomize using a random oracle
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        let m = pk.len();
        let nodes = (0..m).collect::<Vec<_>>();
        let mut positions = vec![];

        for _ in 0..m {
            let mut row = vec![];
            let mut sampled = vec![];
            while sampled.len() < k {
                let value = rng.random_range(0..crs.n);
                // ensure that the sampled value is unique
                if !sampled.contains(&value) {
                    sampled.push(value);
                }
            }
            row.extend(sampled);
            positions.push(row);
        }

        // populate a vector called edges with tuples of (node, position)
        let mut edges = vec![];
        for i in 0..m {
            for j in 0..k {
                edges.push((nodes[i], positions[i][j]));
            }
        }

        use rayon::prelude::*;

        let timer = start_timer!(|| "Setup System Public Keys");
        let mut lag_pk = vec![vec![]; m];
        lag_pk.par_iter_mut().enumerate().for_each(|(i, lag_pk_i)| {
            let mut lag_pk_inner = vec![];
            for j in 0..k {
                lag_pk_inner.push(pk[i].get_lag_public_key(positions[i][j], crs, &lag_polys));
            }
            *lag_pk_i = lag_pk_inner;
            println!("Lagrange public key for party {}", i);
        });
        end_timer!(timer);

        Self { m, k, pk, lag_pk }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::setup::SecretKey;
    use hopcroft_karp::matching;

    type E = ark_bls12_381::Bls12_381;
    type F = ark_bls12_381::Fr;

    #[test]
    fn setup_system_public_keys() {
        let n = 1 << 7;
        let m = 1 << 7;
        let crs = CRS::<E>::new(n, &mut ark_std::test_rng());
        let lag_polys = LagPolys::<F>::new(n);
        use rayon::prelude::*;

        let (_sk, pk): (Vec<_>, Vec<_>) = (0..m)
            .into_par_iter()
            .map(|_| {
                let sk = SecretKey::<E>::new(&mut ark_std::test_rng());
                let pk = sk.get_pk(&crs);
                (sk, pk)
            })
            .unzip();

        let _system_keys = SystemPublicKeys::<E>::new(pk.clone(), &crs, lag_polys, 3);
    }

    #[test]
    fn test_kuhn_munkres() {
        const N: usize = 128;
        const K: usize = 3;
        const RUNS: usize = 1;

        for _ in 0..RUNS {
            let nodes = (0..N).collect::<Vec<_>>();
            // for each node, assign a tuple with K random values from the range 0..128
            let mut rng = rand::rng();
            let mut positions = vec![];
            for _ in 0..N {
                let mut row = vec![];
                let mut sampled = vec![];
                while sampled.len() < K {
                    let value = rng.random_range(N..2 * N);
                    if !sampled.contains(&value) {
                        sampled.push(value);
                    }
                }
                row.extend(sampled);
                positions.push(row);
            }

            // populate a vector called edges with tuples of (node, position)
            let mut edges = vec![];
            for i in 0..N {
                for j in 0..K {
                    edges.push((nodes[i], positions[i][j]));
                }
            }

            let res = matching(&edges);
            println!("Matching Size: {}/{}", res.len(), N);
            println!("Matching: {:?}", res);
        }
    }
}
