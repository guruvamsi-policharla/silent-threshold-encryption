use ark_ec::pairing::Pairing;
use ark_ff::{FftField, Field};
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Evaluations, Polynomial,
    Radix2EvaluationDomain,
};
use ark_std::Zero;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, Validate};

pub(crate) fn ark_se<S, A: CanonicalSerialize>(a: &A, s: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    let mut bytes = vec![];
    a.serialize_with_mode(&mut bytes, Compress::Yes)
        .map_err(serde::ser::Error::custom)?;
    s.serialize_bytes(&bytes)
}

pub(crate) fn ark_de<'de, D, A: CanonicalDeserialize>(data: D) -> Result<A, D::Error>
where
    D: serde::de::Deserializer<'de>,
{
    let s: Vec<u8> = serde::de::Deserialize::deserialize(data)?;
    let a = A::deserialize_with_mode(s.as_slice(), Compress::Yes, Validate::Yes);
    a.map_err(serde::de::Error::custom)
}

// 1 at omega^i and 0 elsewhere on domain {omega^i}_{i \in [n]}
pub fn lagrange_poly<F: FftField>(n: usize, i: usize) -> DensePolynomial<F> {
    debug_assert!(i < n);
    //todo: check n is a power of 2
    let mut evals = vec![];
    for j in 0..n {
        let l_of_x: u64 = if i == j { 1 } else { 0 };
        evals.push(F::from(l_of_x));
    }

    //powers of nth root of unity
    let domain = Radix2EvaluationDomain::<F>::new(n).unwrap();
    let eval_form = Evaluations::from_vec_and_domain(evals, domain);
    //interpolated polynomial over the n points
    eval_form.interpolate()
}

/// interpolates a polynomial which is zero on points and 1 at the point 0
/// todo: use faster interpolation
pub fn interp_mostly_zero<F: Field>(points: &Vec<F>) -> DensePolynomial<F> {
    if points.is_empty() {
        // threshold=n
        return DensePolynomial::from_coefficients_vec(vec![F::one()]);
    }

    let mut interp = DensePolynomial::from_coefficients_vec(vec![F::one()]);
    for &point in points {
        interp = interp.naive_mul(&DensePolynomial::from_coefficients_vec(vec![
            -point,
            F::one(),
        ]));
    }

    let scale = interp.evaluate(&F::zero());
    interp = &interp * (F::one() / scale);

    interp
}

/// Computes all the openings of a KZG commitment in O(n log n) time
/// See https://github.com/khovratovich/Kate/blob/master/Kate_amortized.pdf
/// eprint version has a bug and hasn't been updated
pub fn open_all_values<E: Pairing>(
    y: &[E::G1Affine],
    f: &[E::ScalarField],
    domain: &Radix2EvaluationDomain<E::ScalarField>,
) -> Vec<E::G1> {
    let top_domain = Radix2EvaluationDomain::<E::ScalarField>::new(2 * domain.size()).unwrap();

    // use FK22 to get all the KZG proofs in O(nlog n) time =======================
    // f = {f0 ,f1, ..., fd}
    // v = {(d 0s), f1, ..., fd}
    let mut v = vec![E::ScalarField::zero(); domain.size() + 1];
    v.append(&mut f[1..f.len()].to_vec());

    debug_assert_eq!(v.len(), 2 * domain.size());
    let v = top_domain.fft(&v);

    // h = y \odot v
    let mut h = vec![E::G1::zero(); 2 * domain.size()];
    for i in 0..2 * domain.size() {
        h[i] = y[i] * (v[i]);
    }

    // inverse fft on h
    let mut h = top_domain.ifft(&h);

    h.truncate(domain.size());

    // fft on h to get KZG proofs
    domain.fft(&h)
}

/// interpolates a polynomial where evaluations on points are zero and the polynomial evaluates to 1
/// at the point 1 but relies on the number of points being a power of 2
/// currently not used as this portion is not a bottleneck
// pub fn compute_vanishing_poly(points: &Vec<ScalarField>) -> DensePolynomial {
//     let mut monomials = Vec::new();
//     for i in 0..points.len() {
//         monomials.push(DensePolynomial::from_coeffs(
//             HostSlice::from_slice(&vec![ScalarField::zero() - points[i], ScalarField::one()]),
//             2,
//         ));
//     }

//     // assert that points.len() is a power of 2
//     assert_eq!(
//         points.len().count_ones(),
//         1,
//         "Implementation demands that n-t is a power of 2. Currently: {}",
//         points.len()
//     );

//     let mut chunk_size = points.len() / 2;
//     while chunk_size > 0 {
//         for i in 0..chunk_size {
//             monomials[i] = &monomials[i] * &monomials[i + chunk_size];
//         }
//         chunk_size = chunk_size / 2;
//     }

//     let scale = monomials[0].eval(&ScalarField::one());
//     let res = &monomials[0] * &scale.inv();

//     res
// }

#[cfg(test)]
mod tests {
    use ark_bls12_381::Bls12_381;
    use ark_ec::{bls12::Bls12, pairing::Pairing, VariableBaseMSM};
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
    use ark_std::{UniformRand, Zero};

    use crate::crs::CRS;

    use super::*;
    type Fr = <Bls12<ark_bls12_381::Config> as Pairing>::ScalarField;
    type G1 = <Bls12<ark_bls12_381::Config> as Pairing>::G1;
    type E = Bls12_381;

    #[test]
    fn open_all_test() {
        let mut rng = ark_std::test_rng();

        let n = 1 << 8;
        let domain = Radix2EvaluationDomain::<Fr>::new(n).unwrap();
        let crs = CRS::<E>::new(n, &mut ark_std::test_rng());

        let mut f: Vec<ark_ff::Fp<ark_ff::MontBackend<ark_bls12_381::FrConfig, 4>, 4>> =
            vec![Fr::zero(); n];
        for i in 0..n {
            f[i] = Fr::rand(&mut rng);
        }

        let com = G1::msm(&crs.powers_of_g[0..f.len()], &f).unwrap();

        let timer = std::time::Instant::now();
        let pi = open_all_values::<E>(&crs.y, &f, &domain);
        println!("open_all_values took {:?}", timer.elapsed());

        // verify the kzg proof
        let g = crs.powers_of_g[0];
        let h = crs.powers_of_h[0];

        let fpoly = DensePolynomial::from_coefficients_vec(f.clone());
        for i in 0..n {
            let lhs = E::pairing(com - (g * fpoly.evaluate(&domain.element(i))), h);
            let rhs = E::pairing(pi[i], crs.powers_of_h[1] - (h * domain.element(i)));
            assert_eq!(lhs, rhs);
        }
    }
}
