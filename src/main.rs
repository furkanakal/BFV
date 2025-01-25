use rand::Rng;
use num_bigint::BigInt;
use std::ops::Add;
use num_traits::cast::ToPrimitive;

#[derive(Clone)]
struct BFVParams {
    n: usize,          // Polynomial degree (power of 2)
    t: BigInt,         // Plaintext modulus (prime)
    q: BigInt,         // Ciphertext modulus (prime)
    delta: BigInt,     // Floor(q/t)
} 

#[derive(Clone)]
struct Polynomial {
    coeffs: Vec<BigInt>,
    params: BFVParams,
}

#[derive(Clone)]
struct Ciphertext {
    c0: Polynomial,
    c1: Polynomial,
}

impl BFVParams {
    fn new(n: usize, t: BigInt, q: BigInt) -> Self {
        assert!(n.is_power_of_two(), "n must be a power of 2");
        assert!(t < q, "t must be smaller than q");
        let delta = &q / &t;
        BFVParams { n, t, q, delta }
    }

    fn mod_q(&self, x: &BigInt) -> BigInt {
        let mut v = x % &self.q;
        if v < BigInt::from(0) {
            v += &self.q;
        }
        v
    }
}

impl Polynomial {
    // Binary distribution for secret key
    fn secret_key(params: &BFVParams) -> Self {
        let mut rng = rand::thread_rng();
        let coeffs: Vec<BigInt> = (0..params.n)
            .map(|_| BigInt::from(rng.gen_range(0..=1)))
            .collect();
        Polynomial {
            coeffs,
            params: params.clone(),
        }
    }

    // Very small uniform random for a(X)
    fn uniform_small(params: &BFVParams) -> Self {
        let mut rng = rand::thread_rng();
        let coeffs: Vec<BigInt> = (0..params.n)
            .map(|_| BigInt::from(rng.gen_range(0..=1)))
            .collect();
        Polynomial {
            coeffs,
            params: params.clone(),
        }
    }

    // Tiny error terms
    fn error(params: &BFVParams) -> Self {
        let mut rng = rand::thread_rng();
        let coeffs: Vec<BigInt> = (0..params.n)
            .map(|_| BigInt::from(if rng.gen_bool(0.5) { 0 } else { 1 }))
            .collect();
        Polynomial {
            coeffs,
            params: params.clone(),
        }
    }

    fn add(&self, other: &Polynomial) -> Polynomial {
        let coeffs: Vec<BigInt> = self.coeffs.iter()
            .zip(other.coeffs.iter())
            .map(|(a, b)| self.params.mod_q(&(a + b)))
            .collect();
        Polynomial {
            coeffs,
            params: self.params.clone(),
        }
    }

    fn mul(&self, other: &Polynomial) -> Polynomial {
        let n = self.params.n;
        let mut result = vec![BigInt::from(0); n];
        
        for i in 0..n {
            for j in 0..n {
                let idx = (i + j) % n;
                let sign = if i + j >= n { -1 } else { 1 };
                let term = &self.coeffs[i] * &other.coeffs[j];
                result[idx] = &result[idx] + sign * term;
            }
        }

        let coeffs: Vec<BigInt> = result.iter()
            .map(|x| self.params.mod_q(x))
            .collect();

        Polynomial {
            coeffs,
            params: self.params.clone(),
        }
    }

    fn scalar_mul(&self, scalar: &BigInt) -> Polynomial {
        let coeffs: Vec<BigInt> = self.coeffs.iter()
            .map(|a| self.params.mod_q(&(a * scalar)))
            .collect();
        Polynomial {
            coeffs,
            params: self.params.clone(),
        }
    }
}

impl Add for Ciphertext {
    type Output = Ciphertext;

    fn add(self, other: Ciphertext) -> Ciphertext {
        Ciphertext {
            c0: self.c0.add(&other.c0),
            c1: self.c1.add(&other.c1),
        }
    }
}

impl Ciphertext {
    fn encrypt(params: &BFVParams, message: Vec<i64>, sk: &Polynomial) -> Self {
        // Convert message to polynomial
        let m_coeffs: Vec<BigInt> = message.iter()
            .map(|&x| BigInt::from(x))
            .chain(std::iter::repeat(BigInt::from(0)))
            .take(params.n)
            .collect();

        let m = Polynomial {
            coeffs: m_coeffs,
            params: params.clone(),
        };

        // Sample small polynomials
        let a = Polynomial::uniform_small(params);
        let e = Polynomial::error(params);

        // Compute c0 = a·s + e + Δ·m
        let delta_m = m.scalar_mul(&params.delta);
        let as_poly = a.mul(sk);
        let c0 = as_poly.add(&e).add(&delta_m);
        
        // c1 = -a
        let c1_coeffs = a.coeffs.iter()
            .map(|x| params.mod_q(&(-x)))
            .collect();
        let c1 = Polynomial {
            coeffs: c1_coeffs,
            params: params.clone(),
        };

        Ciphertext { c0, c1 }
    }

    fn decrypt(&self, sk: &Polynomial) -> Vec<i64> {
        // Compute v = c0 + c1·s
        let c1s = self.c1.mul(sk);
        let v = self.c0.add(&c1s);

        // Scale by t/q and round
        let t = &self.c0.params.t;
        let q = &self.c0.params.q;
        
        v.coeffs.iter()
            .map(|x| {
                let scaled = (x * t) / q;
                let val = scaled.to_i64().unwrap_or(0);
                let t_val = t.to_i64().unwrap_or(0);
                let mut result = val % t_val;
                if result < 0 {
                    result += t_val;
                }
                result
            })
            .collect()
    }
}

fn main() {
    let params = BFVParams::new(
        4,                      // n = 4
        BigInt::from(4),       // t = 4 (smaller plaintext modulus)
        BigInt::from(64),      // q = 64 (smaller ciphertext modulus)
    );

    let sk = Polynomial::secret_key(&params);

    // Small plaintext values
    let pt1 = vec![1, 2, 1, 0];
    let pt2 = vec![1, 1, 0, 2];

    println!("Plaintext 1: {:?}", pt1);
    println!("Plaintext 2: {:?}", pt2);

    let ct1 = Ciphertext::encrypt(&params, pt1.clone(), &sk);
    let ct2 = Ciphertext::encrypt(&params, pt2.clone(), &sk);

    println!("\nCiphertext 1:");
    println!("c0: {:?}", ct1.c0.coeffs);
    println!("c1: {:?}", ct1.c1.coeffs);

    println!("\nCiphertext 2:");
    println!("c0: {:?}", ct2.c0.coeffs);
    println!("c1: {:?}", ct2.c1.coeffs);

    let ct_sum = ct1 + ct2;

    println!("\nSum ciphertext:");
    println!("c0: {:?}", ct_sum.c0.coeffs);
    println!("c1: {:?}", ct_sum.c1.coeffs);

    let decrypted = ct_sum.decrypt(&sk);
    println!("\nDecrypted sum: {:?}", decrypted);
    
    let expected: Vec<i64> = pt1.iter()
        .zip(pt2.iter())
        .map(|(&a, &b)| ((a + b) % params.t.to_i64().unwrap()))
        .collect();
    println!("Expected sum:   {:?}", expected);
}