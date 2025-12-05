/// Polynomial operations for threshold encryption
/// Provides functionality similar to ark-poly but using blstrs Scalar field
use blstrs::Scalar;
use ff::{Field, PrimeField};
use std::ops::{Add, Div, Mul, Sub};

/// Dense univariate polynomial represented by coefficients
#[derive(Clone, Debug, PartialEq)]
pub struct DensePolynomial {
    /// Coefficients of the polynomial, where coeffs[i] is the coefficient of x^i
    pub coeffs: Vec<Scalar>,
}

impl DensePolynomial {
    /// Create a new polynomial from coefficients
    pub fn from_coefficients_vec(coeffs: Vec<Scalar>) -> Self {
        let mut poly = DensePolynomial { coeffs };
        poly.truncate_leading_zeros();
        poly
    }

    /// Create zero polynomial
    pub fn zero() -> Self {
        DensePolynomial {
            coeffs: vec![Scalar::ZERO],
        }
    }

    /// Remove leading zero coefficients
    fn truncate_leading_zeros(&mut self) {
        while self.coeffs.len() > 1 && self.coeffs.last() == Some(&Scalar::ZERO) {
            self.coeffs.pop();
        }
        if self.coeffs.is_empty() {
            self.coeffs.push(Scalar::ZERO);
        }
    }

    /// Get the degree of the polynomial
    pub fn degree(&self) -> usize {
        if self.coeffs.len() == 1 && self.coeffs[0] == Scalar::ZERO {
            0
        } else {
            self.coeffs.len().saturating_sub(1)
        }
    }

    /// Evaluate the polynomial at a point using Horner's method
    pub fn evaluate(&self, point: &Scalar) -> Scalar {
        if self.coeffs.is_empty() {
            return Scalar::ZERO;
        }

        let mut result = *self.coeffs.last().unwrap();
        for coeff in self.coeffs.iter().rev().skip(1) {
            result = result * point + coeff;
        }
        result
    }

    /// Naive polynomial multiplication (O(n^2))
    pub fn naive_mul(&self, other: &DensePolynomial) -> DensePolynomial {
        if self.coeffs.is_empty() || other.coeffs.is_empty() {
            return DensePolynomial::zero();
        }

        let mut result = vec![Scalar::ZERO; self.coeffs.len() + other.coeffs.len() - 1];

        for (i, a) in self.coeffs.iter().enumerate() {
            for (j, b) in other.coeffs.iter().enumerate() {
                result[i + j] += *a * b;
            }
        }

        DensePolynomial::from_coefficients_vec(result)
    }

    /// Divide polynomial by vanishing polynomial Z_H(x) = x^n - 1
    /// Returns (quotient, remainder)
    pub fn divide_by_vanishing_poly(
        &self,
        domain: Radix2EvaluationDomain,
    ) -> (DensePolynomial, DensePolynomial) {
        let n = domain.size;

        // Vanishing polynomial is x^n - 1
        if self.degree() < n {
            return (DensePolynomial::zero(), self.clone());
        }

        let mut remainder = self.coeffs.clone();
        let mut quotient = vec![Scalar::ZERO; self.coeffs.len().saturating_sub(n)];

        // Polynomial long division by x^n - 1
        for i in (n..=self.degree()).rev() {
            let coeff = remainder[i];
            quotient[i - n] = coeff;
            // Subtract coeff * (x^n - 1) * x^{i-n} = coeff * x^i - coeff * x^{i-n}
            remainder[i] = Scalar::ZERO;
            remainder[i - n] += coeff;
        }

        remainder.truncate(n);

        let mut quot_poly = DensePolynomial::from_coefficients_vec(quotient);
        let mut rem_poly = DensePolynomial::from_coefficients_vec(remainder);

        quot_poly.truncate_leading_zeros();
        rem_poly.truncate_leading_zeros();

        (quot_poly, rem_poly)
    }
}

impl Add for DensePolynomial {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let max_len = self.coeffs.len().max(other.coeffs.len());
        let mut result = vec![Scalar::ZERO; max_len];

        for (i, coeff) in self.coeffs.iter().enumerate() {
            result[i] = *coeff;
        }
        for (i, coeff) in other.coeffs.iter().enumerate() {
            result[i] += coeff;
        }

        DensePolynomial::from_coefficients_vec(result)
    }
}

impl Sub for DensePolynomial {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let max_len = self.coeffs.len().max(other.coeffs.len());
        let mut result = vec![Scalar::ZERO; max_len];

        for (i, coeff) in self.coeffs.iter().enumerate() {
            result[i] = *coeff;
        }
        for (i, coeff) in other.coeffs.iter().enumerate() {
            result[i] -= coeff;
        }

        DensePolynomial::from_coefficients_vec(result)
    }
}

impl Mul<Scalar> for &DensePolynomial {
    type Output = DensePolynomial;

    fn mul(self, scalar: Scalar) -> DensePolynomial {
        let coeffs = self.coeffs.iter().map(|c| *c * scalar).collect();
        DensePolynomial::from_coefficients_vec(coeffs)
    }
}

impl Mul<&DensePolynomial> for &DensePolynomial {
    type Output = DensePolynomial;

    fn mul(self, other: &DensePolynomial) -> DensePolynomial {
        self.naive_mul(other)
    }
}

impl Div for &DensePolynomial {
    type Output = DensePolynomial;

    fn div(self, divisor: &DensePolynomial) -> DensePolynomial {
        // Polynomial long division
        if divisor.coeffs.len() == 1 && divisor.coeffs[0] == Scalar::ZERO {
            panic!("Division by zero polynomial");
        }

        if self.degree() < divisor.degree() {
            return DensePolynomial::zero();
        }

        let mut remainder = self.clone();
        let mut quotient = vec![Scalar::ZERO; self.degree() - divisor.degree() + 1];

        let divisor_leading_inv = divisor.coeffs.last().unwrap().invert().unwrap();

        for i in (0..=self.degree() - divisor.degree()).rev() {
            let pos = i + divisor.degree();
            let coeff = remainder.coeffs[pos] * divisor_leading_inv;
            quotient[i] = coeff;

            for (j, &div_coeff) in divisor.coeffs.iter().enumerate() {
                remainder.coeffs[i + j] -= coeff * div_coeff;
            }
        }

        DensePolynomial::from_coefficients_vec(quotient)
    }
}

/// Radix-2 FFT evaluation domain for polynomials
#[derive(Clone, Debug)]
pub struct Radix2EvaluationDomain {
    /// Size of the domain (must be a power of 2)
    pub size: usize,
    /// Generator of the domain (n-th root of unity)
    pub group_gen: Scalar,
    /// Inverse of the generator
    pub group_gen_inv: Scalar,
}

impl Radix2EvaluationDomain {
    /// Create a new evaluation domain of the given size
    /// Size must be a power of 2
    pub fn new(size: usize) -> Option<Self> {
        if !size.is_power_of_two() || size == 0 {
            return None;
        }

        // For BLS12-381 scalar field, we need to find the n-th root of unity
        // The scalar field has 2^32 as the two-adicity
        // We compute omega = ROOT_OF_UNITY^{2^32 / size}

        let log_size = size.trailing_zeros();
        if log_size > 32 {
            // BLS12-381 scalar field has 2-adicity of 32
            return None;
        }

        // Get the primitive 2^32-th root of unity for BLS12-381
        // This is a constant for the BLS12-381 scalar field
        let root_of_unity = Scalar::ROOT_OF_UNITY;

        // Compute the n-th root of unity: omega = root^{2^32 / n}
        // We need to compute root_of_unity^{2^{32-log_size}}
        let power = 1u64 << (32 - log_size);
        let mut group_gen = Scalar::ONE;
        let mut base = root_of_unity;
        let mut exp = power;

        // Fast exponentiation
        while exp > 0 {
            if exp & 1 == 1 {
                group_gen *= base;
            }
            base = base * base;
            exp >>= 1;
        }

        let group_gen_inv = group_gen.invert().unwrap();

        Some(Radix2EvaluationDomain {
            size,
            group_gen,
            group_gen_inv,
        })
    }

    /// Return an iterator over the elements of the domain
    pub fn elements(&self) -> impl Iterator<Item = Scalar> + '_ {
        let mut current = Scalar::ONE;
        let generator = self.group_gen;
        (0..self.size).map(move |i| {
            let result = current;
            if i < self.size - 1 {
                current *= generator;
            }
            result
        })
    }

    /// Perform FFT (Fast Fourier Transform) on polynomial coefficients
    /// Converts from coefficient representation to evaluation representation
    pub fn fft(&self, coeffs: &[Scalar]) -> Vec<Scalar> {
        let mut a = coeffs.to_vec();
        a.resize(self.size, Scalar::ZERO);

        self.fft_in_place(&mut a);
        a
    }

    /// In-place FFT implementation
    fn fft_in_place(&self, a: &mut [Scalar]) {
        let n = a.len();
        assert!(n == self.size);

        if n == 1 {
            return;
        }

        // Bit-reversal permutation
        let mut j = 0;
        for i in 1..n {
            let mut bit = n >> 1;
            while j & bit != 0 {
                j ^= bit;
                bit >>= 1;
            }
            j ^= bit;

            if i < j {
                a.swap(i, j);
            }
        }

        // Cooley-Tukey FFT
        let mut len = 2;
        while len <= n {
            let half_len = len / 2;

            // Compute omega for this stage
            let angle = self.size / len;
            let mut omega_step = Scalar::ONE;
            for _ in 0..angle {
                omega_step *= self.group_gen;
            }

            for start in (0..n).step_by(len) {
                let mut omega = Scalar::ONE;
                for j in 0..half_len {
                    let u = a[start + j];
                    let v = a[start + j + half_len] * omega;
                    a[start + j] = u + v;
                    a[start + j + half_len] = u - v;
                    omega *= omega_step;
                }
            }

            len *= 2;
        }
    }

    /// Perform inverse FFT
    /// Converts from evaluation representation to coefficient representation
    pub fn ifft(&self, evals: &[Scalar]) -> Vec<Scalar> {
        let mut a = evals.to_vec();
        a.resize(self.size, Scalar::ZERO);

        self.ifft_in_place(&mut a);
        a
    }

    /// In-place inverse FFT
    fn ifft_in_place(&self, a: &mut [Scalar]) {
        let n = a.len();

        // Use the inverse generator
        let mut domain_inv = self.clone();
        domain_inv.group_gen = self.group_gen_inv;

        domain_inv.fft_in_place(a);

        // Scale by 1/n
        let n_inv = Scalar::from(n as u64).invert().unwrap();
        for coeff in a.iter_mut() {
            *coeff *= n_inv;
        }
    }
}

/// Polynomial in evaluation form
#[derive(Clone, Debug)]
pub struct Evaluations {
    pub evals: Vec<Scalar>,
    pub domain: Radix2EvaluationDomain,
}

impl Evaluations {
    pub fn from_vec_and_domain(evals: Vec<Scalar>, domain: Radix2EvaluationDomain) -> Self {
        assert_eq!(evals.len(), domain.size);
        Evaluations { evals, domain }
    }

    /// Interpolate the polynomial from its evaluations using inverse FFT
    pub fn interpolate(self) -> DensePolynomial {
        let coeffs = self.domain.ifft(&self.evals);
        DensePolynomial::from_coefficients_vec(coeffs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polynomial_evaluation() {
        // p(x) = 1 + 2x + 3x^2
        let poly = DensePolynomial::from_coefficients_vec(vec![
            Scalar::ONE,
            Scalar::from(2),
            Scalar::from(3),
        ]);

        // p(0) = 1
        assert_eq!(poly.evaluate(&Scalar::ZERO), Scalar::ONE);

        // p(1) = 1 + 2 + 3 = 6
        assert_eq!(poly.evaluate(&Scalar::ONE), Scalar::from(6));
    }

    #[test]
    fn test_polynomial_addition() {
        let p1 = DensePolynomial::from_coefficients_vec(vec![Scalar::ONE, Scalar::from(2)]);
        let p2 = DensePolynomial::from_coefficients_vec(vec![Scalar::from(3), Scalar::from(4)]);
        let sum = p1 + p2;

        assert_eq!(sum.coeffs, vec![Scalar::from(4), Scalar::from(6)]);
    }

    #[test]
    fn test_polynomial_multiplication() {
        // (1 + x) * (1 + x) = 1 + 2x + x^2
        let p = DensePolynomial::from_coefficients_vec(vec![Scalar::ONE, Scalar::ONE]);
        let result = p.naive_mul(&p);

        assert_eq!(result.coeffs.len(), 3);
        assert_eq!(result.coeffs[0], Scalar::ONE);
        assert_eq!(result.coeffs[1], Scalar::from(2));
        assert_eq!(result.coeffs[2], Scalar::ONE);
    }

    #[test]
    fn test_fft_domain() {
        let domain = Radix2EvaluationDomain::new(4).unwrap();
        let elements: Vec<_> = domain.elements().collect();

        assert_eq!(elements.len(), 4);
        // omega^4 should equal 1
        let omega4 = elements[1] * elements[1] * elements[1] * elements[1];
        assert_eq!(omega4, Scalar::ONE);
    }
}
