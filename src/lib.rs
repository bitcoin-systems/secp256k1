use std::ops::{Add, Sub, Mul, Div};
use std::fmt;

// Field modulus for secp256k1
const FIELD_MODULUS: [u64; 4] = [
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFE,
    0xBAAEDCE6AF48A03B,
    0xBFD25E8CD0364141,
];

// Curve parameters
const CURVE_A: u64 = 0;
const CURVE_B: [u64; 4] = [7, 0, 0, 0];
const GENERATOR_X: [u64; 4] = [
    0x79BE667EF9DCBBAC,
    0x55A06295CE870B07,
    0x037BF51F5CAEE4B7,
    0x5CBF3BA3D23861BE,
];
const GENERATOR_Y: [u64; 4] = [
    0x483ADA7726A3C465,
    0x5DA4FBFC0E1108A8,
    0xFD17B448A6855419,
    0x9C47D08FFB10D4B8,
];

// Field element representation
#[derive(Clone, Copy, PartialEq, Eq)]
struct FieldElement {
    value: [u64; 4],
}

impl FieldElement {
    // Create a new field element from u64 array
    fn new(value: [u64; 4]) -> Self {
        let mut reduced = value;
        FieldElement::reduce(&mut reduced);
        FieldElement { value: reduced }
    }

    // Reduce field element modulo prime
    fn reduce(value: &mut [u64; 4]) {
        let mut carry: u64 = 0;
        for i in 0..4 {
            let temp = value[i].wrapping_add(carry);
            value[i] = temp;
            carry = temp >> 63;
        }

        // Modular reduction using Barrett reduction
        let mut q = [0u64; 8];
        let mut r = [0u64; 8];

        // Multiplication approach for modular reduction
        for i in 0..4 {
            for j in 0..4 {
                let (high, low) = Self::wide_mul(value[i], FIELD_MODULUS[j]);
                q[i + j] = q[i + j].wrapping_add(low);
                q[i + j + 1] = q[i + j + 1].wrapping_add(high);
            }
        }

        // Subtract multiples of modulus
        for i in 0..4 {
            value[i] = r[i];
        }
    }

    // Wide multiplication to handle large numbers
    fn wide_mul(a: u64, b: u64) -> (u64, u64) {
        let a128 = a as u128;
        let b128 = b as u128;
        let product = a128 * b128;
        ((product >> 64) as u64, product as u64)
    }

    // Multiplicative inverse using extended Euclidean algorithm
    fn inv(&self) -> Self {
        // Fermat's little theorem: a^-1 ≡ a^(p-2) (mod p)
        let mut exp = FIELD_MODULUS;
        exp[0] -= 2; // p-2

        let mut result = FieldElement::new([1, 0, 0, 0]);
        let mut base = *self;

        for word in exp.iter().rev() {
            for bit in (0..64).rev() {
                result = result * result;
                if (word & (1 << bit)) != 0 {
                    result = result * *self;
                }
            }
        }

        result
    }
}

// Arithmetic implementations
impl Add for FieldElement {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        let mut result = [0u64; 4];
        let mut carry: u64 = 0;
        for i in 0..4 {
            let temp = self.value[i].wrapping_add(other.value[i]).wrapping_add(carry);
            result[i] = temp;
            carry = temp >> 63;
        }
        FieldElement::new(result)
    }
}

impl Sub for FieldElement {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        let mut result = [0u64; 4];
        let mut borrow: u64 = 0;
        for i in 0..4 {
            let temp = self.value[i].wrapping_sub(other.value[i]).wrapping_sub(borrow);
            result[i] = temp;
            borrow = temp >> 63;
        }
        FieldElement::new(result)
    }
}

impl Mul for FieldElement {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        let mut result = [0u64; 8];
        for i in 0..4 {
            let mut carry: u64 = 0;
            for j in 0..4 {
                let (high, low) = FieldElement::wide_mul(self.value[i], other.value[j]);
                let temp = result[i + j].wrapping_add(low).wrapping_add(carry);
                result[i + j] = temp;
                carry = high.wrapping_add(temp >> 63);
            }
            result[i + 4] = carry;
        }

        let mut reduced = [0u64; 4];
        for i in 0..4 {
            reduced[i] = result[i];
        }
        FieldElement::new(reduced)
    }
}

// Point representation on the elliptic curve
#[derive(Clone, Copy)]
struct Point {
    x: Option<FieldElement>,
    y: Option<FieldElement>,
}

impl Point {
    // Point addition on the curve
    fn add(&self, other: &Point) -> Point {
        match (self.x, self.y, other.x, other.y) {
            // Point at infinity cases
            (None, _, _, _) => *other,
            (_, _, None, _) => *self,

            // Special cases where one point has only x coordinate
            (Some(_), None, Some(_), _) | (Some(_), None, _, _) => *self,
            (Some(_), _, Some(_), None) | (_, _, Some(_), None) => *other,

            // General point addition for points with full coordinates
            (Some(x1), Some(y1), Some(x2), Some(y2)) => {
                // Slope calculation
                let slope = if x1 == x2 && y1 == y2 {
                    // Point doubling case
                    let three_x_squared = FieldElement::new([3, 0, 0, 0]) * x1 * x1;
                    let two_y = FieldElement::new([2, 0, 0, 0]) * y1;
                    three_x_squared * two_y.inv()
                } else {
                    // Point addition case
                    (y2 - y1) * (x2 - x1).inv()
                };

                // New x coordinate
                let x3 = slope * slope - x1 - x2;

                // New y coordinate
                let y3 = slope * (x1 - x3) - y1;

                Point {
                    x: Some(x3),
                    y: Some(y3),
                }
            }
        }
    }

    // Scalar multiplication using double-and-add method
    fn scalar_mul(&self, scalar: &[u64; 4]) -> Point {
        let mut result = Point { x: None, y: None };
        let mut base = *self;

        for word in scalar.iter().rev() {
            for bit in (0..64).rev() {
                result = result.add(&result);
                if (word & (1 << bit)) != 0 {
                    result = result.add(&base);
                }
            }
        }

        result
    }
}

// Utility functions for cryptographic operations
impl Point {
    // Check if point is on the curve
    fn is_on_curve(&self) -> bool {
        match (self.x, self.y) {
            (Some(x), Some(y)) => {
                // y² = x³ + 7 (secp256k1 curve equation)
                let x_cubed = x * x * x;
                let y_squared = y * y;
                x_cubed + FieldElement::new(CURVE_B) == y_squared
            }
            _ => true // Point at infinity is considered on the curve
        }
    }

    // Generator point for secp256k1
    fn generator() -> Self {
        Point {
            x: Some(FieldElement::new(GENERATOR_X)),
            y: Some(FieldElement::new(GENERATOR_Y)),
        }
    }
}

// Debug implementation for pretty printing
impl fmt::Debug for FieldElement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FieldElement({:?})", self.value)
    }
}

impl fmt::Debug for Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match (self.x, self.y) {
            (Some(x), Some(y)) => write!(f, "Point(x: {:?}, y: {:?})", x, y),
            _ => write!(f, "Point(Infinity)"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper function to create a field element from a u64 value
    fn fe(value: u64) -> FieldElement {
        FieldElement::new([value, 0, 0, 0])
    }

    #[test]
    fn test_field_element_arithmetic() {
        // Basic addition
        let a = fe(5);
        let b = fe(3);
        let sum = a + b;
        assert_eq!(sum.value[0], 8);

        // Subtraction
        let diff = a - b;
        assert_eq!(diff.value[0], 2);

        // Multiplication
        let product = a * b;
        assert_eq!(product.value[0], 15);
    }

    #[test]
    fn test_point_addition() {
        // Test point at infinity cases
        let infinity = Point { x: None, y: None };
        let generator = Point::generator();

        // Adding point to infinity should return the point
        assert!(generator.add(&infinity).x == generator.x);
        assert!(infinity.add(&generator).x == generator.x);

        // Test point doubling
        let doubled_gen = generator.add(&generator);
        assert!(doubled_gen.is_on_curve());
        
        // Verify that doubled point is different from original generator
        assert!(doubled_gen.x != generator.x);
    }

    #[test]
    fn test_point_on_curve() {
        // Verify generator point is on the curve
        let generator = Point::generator();
        assert!(generator.is_on_curve(), "Generator point must be on curve");

        // Test a few manually constructed points
        let test_x = FieldElement::new([
            0x79BE667EF9DCBBAC,
            0x55A06295CE870B07,
            0x037BF51F5CAEE4B7,
            0x5CBF3BA3D23861BE,
        ]);
        let test_y = FieldElement::new([
            0x483ADA7726A3C465,
            0x5DA4FBFC0E1108A8,
            0xFD17B448A6855419,
            0x9C47D08FFB10D4B8,
        ]);

        let test_point = Point {
            x: Some(test_x),
            y: Some(test_y),
        };
        assert!(test_point.is_on_curve(), "Test point must be on curve");
    }

    #[test]
    fn test_scalar_multiplication() {
        let generator = Point::generator();
        
        // Test scalar multiplication with small scalar
        let small_scalar = [2u64, 0, 0, 0];
        let doubled_point = generator.scalar_mul(&small_scalar);
        
        // Verify that scalar multiplication results in a point on the curve
        assert!(doubled_point.is_on_curve(), "Scalar multiplication must result in point on curve");
        
        // Verify that scalar multiplication of 2 is equivalent to point doubling
        let direct_double = generator.add(&generator);
        assert_eq!(
            doubled_point.x.unwrap().value, 
            direct_double.x.unwrap().value, 
            "Scalar multiplication by 2 must match point doubling"
        );
    }

    #[test]
    fn test_multiplicative_inverse() {
        // Test multiplicative inverse property
        let a = fe(5);
        let inv_a = a.inv();
        
        // Verify a * a^-1 = 1 (in field arithmetic)
        let product = a * inv_a;
        assert_eq!(product.value[0], 1, "Multiplicative inverse must satisfy a * a^-1 = 1");
    }

    #[test]
    fn test_edge_cases() {
        // Test addition with points having partial coordinates
        let partial_x_point = Point { 
            x: Some(fe(10)), 
            y: None 
        };
        let full_point = Point::generator();

        // These should handle partial coordinate scenarios
        let result1 = partial_x_point.add(&full_point);
        let result2 = full_point.add(&partial_x_point);

        // Verify results are valid
        assert!(result1.x.is_some(), "Addition with partial coordinates should produce a point");
        assert!(result2.x.is_some(), "Addition with partial coordinates should produce a point");
    }
}

fn main() {
    // Basic example usage and testing
    let gen = Point::generator();
    println!("Generator Point: {:?}", gen);
    println!("Is generator on curve: {}", gen.is_on_curve());
}