//! Common small functions used throughout the crate
//!
//! These are left public for the convenience of the user. For example
//! capitalising a string or using prettier formatting for scientific numbers.

use std::fmt::LowerExp;

// Alias for the format! macro out of laziness
pub use std::format as f;

/// Extends primitives with more specific formatting options
pub trait NumberFmt {
    /// Better scientific number formatting
    ///
    /// The default is not very consistent for scientific in particular, so this
    /// allows easy definition.
    ///
    /// Works for anything that can be represented as scientific using the  
    /// LowerExp trait.
    ///
    /// ```rust
    /// # use meshtal::utils::NumberFmt;
    /// let number = -1.0;
    /// assert_eq!(number.sci(5, 2), "-1.00000e+00".to_string());
    /// assert_eq!((1.0).sci(5, 2), "1.00000e+00".to_string());
    /// ```
    fn sci(&self, precision: usize, exp_pad: usize) -> String;
}

impl<T: LowerExp> NumberFmt for T {
    fn sci(&self, precision: usize, exp_pad: usize) -> String {
        let mut num = f!("{:.precision$e}", &self, precision = precision);
        // Safe to `unwrap` as `num` is guaranteed to contain `'e'`
        let exp = num.split_off(num.find('e').unwrap());
        // Make sure the exponent is signed
        let (sign, exp) = match exp.strip_prefix("e-") {
            Some(exp) => ('-', exp),
            None => ('+', &exp[1..]),
        };
        // Pad the exponent with zeros if needed and put it back on the number
        num.push_str(&f!("e{}{:0>pad$}", sign, exp, pad = exp_pad));
        num
    }
}

/// Capilalises the first letter in a string
///
/// ```rust
/// # use meshtal::utils::capitalise;
/// assert_eq!(capitalise("test string"), "Test string".to_string());
/// ```
pub fn capitalise(s: &str) -> String {
    let mut c = s.chars();
    match c.next() {
        Some(f) => f.to_uppercase().collect::<String>() + c.as_str(),
        None => String::new(),
    }
}

/// Find the maximum value of a `Vec<f64>`
///
/// Rust only havs a built-in max method for types that implement Ord. However,
/// floating-point types do not implement Ord because of NaN, so this is the
/// workaround.
///
/// ```rust
/// # use meshtal::utils::vec_f64_max;
/// let vector = vec![1.0, 2.0, 3.0];
/// assert_eq!(*vec_f64_max(&vector), 3.0)
/// ```
pub fn vec_f64_max(vector: &[f64]) -> &f64 {
    vector.iter().max_by(|a, b| a.total_cmp(b)).unwrap()
}

/// Find the minimum value of a `Vec<f64>`
///
/// Rust only havs a built-in max method for types that implement Ord. However,
/// floating-point types do not implement Ord because of NaN, so this is the
/// workaround.
///
/// ```rust
/// # use meshtal::utils::vec_f64_min;
/// let vector = vec![1.0, 2.0, 3.0];
/// assert_eq!(*vec_f64_min(&vector), 1.0)
/// ```
pub fn vec_f64_min(vector: &[f64]) -> &f64 {
    vector.iter().min_by(|a, b| a.total_cmp(b)).unwrap()
}
