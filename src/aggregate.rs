use crate::{
    crs::CRS,
    setup::{LagPolys, LagPublicKey, PublicKey},
    utils::{ark_de, ark_se},
};
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{end_timer, start_timer, One, Zero};
use hopcroft_karp::matching;
use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

#[derive(CanonicalSerialize, CanonicalDeserialize, Serialize, Deserialize, Clone)]
pub struct EncryptionKey<E: Pairing> {
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub ask: E::G1,
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub z_g2: E::G2,
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub e_gh: PairingOutput<E>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Serialize, Deserialize, Clone)]
pub struct AggregateKey<E: Pairing> {
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub lag_pks: Vec<LagPublicKey<E>>,
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub agg_sk_li_lj_z: Vec<E::G1>,
}

impl<E: Pairing> AggregateKey<E> {
    pub fn new(lag_pks: Vec<LagPublicKey<E>>, crs: &CRS<E>) -> (Self, EncryptionKey<E>) {
        let n = lag_pks.len();
        let z_g2 = crs.powers_of_h[n] + crs.powers_of_h[0] * (-E::ScalarField::one());

        // gather sk_li from all public keys
        let mut ask = E::G1::zero();
        for pki in lag_pks.iter() {
            ask += pki.sk_li;
        }

        let mut agg_sk_li_lj_z = vec![];
        for i in 0..n {
            let mut agg_sk_li_lj_zi = E::G1::zero();
            for pkj in lag_pks.iter() {
                agg_sk_li_lj_zi += pkj.sk_li_lj_z[i];
            }
            agg_sk_li_lj_z.push(agg_sk_li_lj_zi);
        }

        (
            AggregateKey {
                lag_pks,
                agg_sk_li_lj_z,
            },
            EncryptionKey {
                ask,
                z_g2,
                e_gh: E::pairing(crs.powers_of_g[0], crs.powers_of_h[0]),
            },
        )
    }
}

/// contains public keys of all m parties in the system
/// also maintains lagrange public keys for each party at k "random" positions
/// this ensures that the public keys of any subset of n parties forms an "almost" perfect matching
/// hence, allowing for efficient aggregation and decryption
#[derive(CanonicalSerialize, CanonicalDeserialize, Serialize, Deserialize, Clone)]
pub struct SystemPublicKeys<E: Pairing> {
    pub m: usize,
    pub k: usize,
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub pks: Vec<PublicKey<E>>,
    #[serde(serialize_with = "ark_se", deserialize_with = "ark_de")]
    pub lag_pks: Vec<Vec<LagPublicKey<E>>>,
}

impl<E: Pairing> SystemPublicKeys<E> {
    pub fn new(
        pks: Vec<PublicKey<E>>,
        crs: &CRS<E>,
        lag_polys: &LagPolys<E::ScalarField>,
        k: usize,
    ) -> Self {
        // using a deterministic seed for reproducibility across machines
        // can derandomize using a random oracle
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        let m = pks.len();
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
        let mut lag_pks = vec![vec![]; m];
        lag_pks
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, lag_pk_i)| {
                let mut lag_pk_inner = vec![];
                for j in 0..k {
                    lag_pk_inner.push(pks[i].get_lag_public_key(positions[i][j], crs, &lag_polys));
                }
                *lag_pk_i = lag_pk_inner;
            });
        end_timer!(timer);

        Self { m, k, pks, lag_pks }
    }

    pub fn get_aggregate_key(
        &self,
        set: &Vec<usize>, // subset of indexes to encrypt to
        crs: &CRS<E>,
        lag_polys: &LagPolys<E::ScalarField>,
    ) -> (AggregateKey<E>, EncryptionKey<E>) {
        // find the maximum matching between lag_public keys and positions
        let mut edges = vec![];
        for &node in set {
            for j in 0..self.k {
                edges.push((node, self.m + self.lag_pks[node][j].position));
            }
        }

        let mut res = matching(&edges);

        println!("Initial Matching Size: {}/{}", res.len(), crs.n);

        // deterministically assign the remaining edges
        let mut assigned_nodes = vec![];
        let mut assigned_positions = vec![];

        for &(node, position) in res.iter() {
            assigned_nodes.push(node);
            assigned_positions.push(position);
        }

        assigned_nodes.sort();
        assigned_positions.sort();

        let mut unassigned_positions: Vec<_> = (self.m..self.m + crs.n)
            .filter(|position| assigned_positions.binary_search(position).is_err())
            .collect();

        // iterate over the set an assign any unassigned nodes
        for &node in set {
            if assigned_nodes.binary_search(&node).is_err() {
                res.push((node, unassigned_positions[0]));
                unassigned_positions.remove(0);
            }
        }

        // sort res based on the second element of the tuple
        res.sort_by_key(|&(_, position)| position);
        // subtract m from the second element of the tuple
        res.iter_mut().for_each(|&mut (_, ref mut position)| {
            *position -= self.m;
        });

        println!("Matching Size: {}/{}", res.len(), crs.n);
        // create a new vector of lag public keys
        let mut set_lag_pks = vec![];
        for i in 0..res.len() {
            let (node, position) = res[i];
            // check if the position is already present in the lag_pk
            // and if so, push that lag_pk to the agg_pk
            // otherwise, create a new lag public key
            // println!("Node: {}, Position: {}", node, position);
            if let Some(lag_pk) = self.lag_pks[node]
                .iter()
                .find(|&lag_pk| lag_pk.position == position)
            {
                set_lag_pks.push(lag_pk.clone());
            } else {
                set_lag_pks.push(self.pks[node].get_lag_public_key(position, crs, &lag_polys));
            }
        }

        AggregateKey::new(set_lag_pks, crs)
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
        let m = 1 << 10;
        let crs = CRS::<E>::new(n, &mut ark_std::test_rng());
        let lag_polys = LagPolys::<F>::new(n);
        use rayon::prelude::*;

        let timer = start_timer!(|| "Setup Public Keys");
        let (_sk, pk): (Vec<_>, Vec<_>) = (0..m)
            .into_par_iter()
            .map(|i| {
                let sk = SecretKey::<E>::new(&mut ark_std::test_rng(), i);
                let pk = sk.get_pk(&crs);
                (sk, pk)
            })
            .unzip();
        end_timer!(timer);

        let _system_keys = SystemPublicKeys::<E>::new(pk.clone(), &crs, &lag_polys, 3);
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

            let mut res = matching(&edges);

            println!("Initial Matching Size: {}/{}", res.len(), N);

            let mut assigned_nodes = vec![];
            let mut assigned_positions = vec![];

            for &(node, position) in res.iter() {
                assigned_nodes.push(node);
                assigned_positions.push(position);
            }

            assigned_nodes.sort();
            assigned_positions.sort();

            let mut unassigned_positions: Vec<_> = (N..2 * N)
                .filter(|position| assigned_positions.binary_search(position).is_err())
                .collect();

            // iterate over the set an assign any unassigned nodes
            for i in 0..N {
                if assigned_nodes.binary_search(&i).is_err() {
                    res.push((i, unassigned_positions[0]));
                    unassigned_positions.remove(0);
                }
            }

            println!("Matching Size: {}/{}", res.len(), N);
            println!("Matching: {:?}", res);
        }
    }
}
