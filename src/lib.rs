//! # logsumexp
//!
//! Numerically stable evaluation of `log(exp(a) + exp(b))` via the `LogAddExp` trait,
//! and a numerically stable, 1-pass algorithm for evaluation of [LogSumExp](https://en.wikipedia.org/wiki/LogSumExp)
//! via the `LogSumExp` trait.

use lnexp::LnExp;

/// A trait which, for the type on which it is implemented,
/// provides numerically-stable evaluation of `log(exp(a) + exp(b))`.
/// The provided implementations on `f64` and `f32` utilize [`ln_1p_exp`](https://docs.rs/lnexp/0.2.0/lnexp/trait.LnExp.html#tymethod.ln_1p_exp)
/// for maximum stability.
pub trait LogAddExp<Rhs = Self> {
    type Output;

    /// Return the log of the sum of exponentials of `self` and `rhs`.
    ///
    /// # Examples
    /// ```
    /// use logsumexp::LogAddExp;
    ///
    /// let x: f64 = 0.5;
    /// let y: f64 = 1.0;
    /// let z: f64 = (x.exp() + y.exp()).ln();
    /// assert_eq!(x.ln_add_exp(y), z);
    /// assert_eq!(x.ln_add_exp(&y), z);
    /// assert_eq!((&x).ln_add_exp(y), z);
    /// assert_eq!((&x).ln_add_exp(&y), z);
    ///
    /// let x: f64 = 1023.0;
    /// let y: f64 = 511.0;
    /// assert_eq!(x.ln_add_exp(y), x);
    /// // compare to naive computation
    /// assert_eq!((x.exp() + y.exp()).ln(), f64::INFINITY);
    /// ```
    fn ln_add_exp(&self, rhs: Rhs) -> Self::Output;
}

macro_rules! impl_logaddexp {
    { $($f:ident)+ } => {
        $(
            impl LogAddExp for $f {
                type Output = $f;
                fn ln_add_exp(&self, rhs: Self) -> Self::Output {
                    let (max, diff) = if *self < rhs {
                        (rhs, *self - rhs)
                    } else {
                        if *self == rhs {
                            (rhs, 0.0)
                        } else {
                            (*self, rhs - *self)
                        }
                    };
                    max + diff.ln_1p_exp()
                }
            }
            impl LogAddExp<&$f> for $f {
                type Output = $f;
                fn ln_add_exp(&self, rhs: &$f) -> Self::Output {
                    self.ln_add_exp(*rhs)
                }
            }
        )+

    };
}
impl_logaddexp! { f64 f32 }

/// A trait for computing the log of the sum of exponentials of a sequence
/// in a numerically-stable manner, using a 1-pass algorithm based on
/// [Milakov, Maxim, and Natalia Gimelshein. "Online normalizer calculation for softmax." arXiv preprint arXiv:1805.02867 (2018)](https://arxiv.org/pdf/1805.02867.pdf).
/// In contrast to the original, this algorithm correctly handles +/-infinity and `nan` values
/// at any point in the sequence.
pub trait LogSumExp<T, U: Iterator<Item = T>> {
    type Output;

    /// Return the [LogSumExp](https://en.wikipedia.org/wiki/LogSumExp) of the sequence.
    ///
    /// # Examples
    /// ```
    /// use logsumexp::LogSumExp;
    ///
    /// let v: Vec<f64> = [0.5_f64, 0.3, 0.1, 0.7].into_iter().map(|x| x.ln()).collect();
    /// let rhs: f64 = (1.6_f64).ln();
    /// assert!((v.iter().ln_sum_exp() -  rhs).abs() < f64::EPSILON);
    ///
    /// let v = (0..10).into_iter().map(|x| (x as f64).ln());
    /// let rhs: f64 = (45.0_f64).ln();
    /// assert_eq!(v.ln_sum_exp(), rhs);
    ///
    /// // mean on the log scale; useful when working with many small (log-)probabilities
    /// let v: Vec<f64> = vec![0.1, 0.2, 0.3, 0.4, 0.5];
    /// let n = v.len();
    /// let rhs: f64 = 0.3;
    /// let log_mean = v.into_iter().map(|x| x.ln()).ln_sum_exp() - (n as f64).ln();
    /// assert!((log_mean.exp() - rhs).abs() < f64::EPSILON);
    /// ```
    fn ln_sum_exp(self) -> Self::Output;
}

macro_rules! impl_logsumexp {
    { $($f:ident)+ } => {
        $(
            impl<U> LogSumExp<$f, U> for U
            where
                U: Iterator<Item = $f>,
            {
                type Output = $f;
                fn ln_sum_exp(mut self) -> Self::Output {
                    let mut m_old = $f::NEG_INFINITY;
                    let mut sum: $f = 0.0;
                    while let Some(v_i) = self.next() {
                        // This is the concept, but it can probably invoke fewer branches.
                        if v_i == $f::NEG_INFINITY {
                            // Of the special cases, -inf is the most likely, hence,
                            // check for it first.
                            continue
                        } else if v_i == $f::INFINITY {
                            // inf should be more likely than nan, under reasonable
                            // circumstances.
                            while let Some(v_i) = self.next() {
                                if v_i.is_nan() {
                                    return v_i
                                }
                            }
                            return $f::INFINITY
                        } else if v_i.is_nan() {
                            // The check for nan is unavoidable.
                            return v_i
                        } else {
                            // finite and not nan
                            let m_new = m_old.max(v_i);
                            sum = sum * (m_old - m_new).exp() + (v_i - m_new).exp();
                            m_old = m_new;
                        }
                    }
                    m_old + sum.ln()
                }
            }

            impl<'a, U> LogSumExp<&'a $f, U> for U
            where
                U: Iterator<Item = &'a $f>,
            {
                type Output = $f;
                fn ln_sum_exp(mut self) -> Self::Output {
                    let mut m_old = $f::NEG_INFINITY;
                    let mut sum: $f = 0.0;
                    while let Some(v_i) = self.next() {
                        if *v_i == $f::NEG_INFINITY {
                            continue
                        } else if *v_i == $f::INFINITY {
                            while let Some(v_i) = self.next() {
                                if v_i.is_nan() {
                                    return *v_i
                                }
                            }
                            return $f::INFINITY
                        } else if v_i.is_nan() {
                            return *v_i
                        } else {
                            let m_new = m_old.max(*v_i);
                            sum = sum * (m_old - m_new).exp() + (*v_i - m_new).exp();
                            m_old = m_new;
                        }
                    }
                    m_old + sum.ln()
                }
            }
        )+

    }
}
impl_logsumexp! { f64 f32 }

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! ln_add_exp_tests {
        { $name:ident $f:ident } => {
            #[cfg(test)]
            mod $name {
                use super::*;
                #[test]
                fn ln_add_exp_works() {
                    let inf = $f::INFINITY;
                    let nan: $f = $f::NAN;

                    // Cases involving +/-inf
                    let u: $f = inf;
                    let v: $f = inf;
                    assert_eq!(u.ln_add_exp(v), inf);

                    let u: $f = -inf;
                    let v: $f = -inf;
                    assert_eq!(u.ln_add_exp(v), -inf);

                    let u: $f = inf;
                    let v: $f = -inf;
                    assert_eq!(u.ln_add_exp(v), inf);

                    let u: $f = -inf;
                    let v: $f = inf;
                    assert_eq!(u.ln_add_exp(v), inf);

                    let u: $f = inf;
                    let v: $f = 0.5;
                    assert_eq!(u.ln_add_exp(v), inf);

                    let u: $f = -inf;
                    let v: $f = 0.5;
                    assert_eq!(u.ln_add_exp(v), v);

                    let u: $f = 0.5;
                    let v: $f = inf;
                    assert_eq!(u.ln_add_exp(v), inf);

                    let u: $f = 0.5;
                    let v: $f = -inf;
                    assert_eq!(u.ln_add_exp(v), u);

                    // Cases involving nan
                    assert!(nan.ln_add_exp(inf).is_nan());
                    assert!(nan.ln_add_exp(-inf).is_nan());
                    assert!(inf.ln_add_exp(nan).is_nan());
                    assert!((-inf).ln_add_exp(nan).is_nan());

                    let x: $f = 0.5;
                    assert!(nan.ln_add_exp(x).is_nan());
                    assert!(nan.ln_add_exp(-x).is_nan());
                    assert!(x.ln_add_exp(nan).is_nan());
                    assert!((-x).ln_add_exp(nan).is_nan());
                }

                #[test]
                fn ln_add_exp_works_argtypes() {
                    let x: $f = 0.5;
                    let y: $f = 1.0;
                    let z: $f = (x.exp() + y.exp()).ln();
                    assert_eq!(x.ln_add_exp(y), z);
                    assert_eq!(x.ln_add_exp(&y), z);
                    let x_ref = &x;
                    assert_eq!(x_ref.ln_add_exp(y), z);
                    assert_eq!(x_ref.ln_add_exp(&y), z);
                }
            }
        }
    }

    ln_add_exp_tests! { f64_logaddexp_impl f64 }
    ln_add_exp_tests! { f32_logaddexp_impl f32 }

    macro_rules! ln_sum_exp_tests {
        { $name:ident $f:ident } => {
            #[cfg(test)]
            mod $name {
                use super::*;

                #[test]
                fn ln_sum_exp_works() {
                    let inf: $f = $f::INFINITY;
                    let neg_inf: $f = $f::NEG_INFINITY;
                    let nan: $f = $f::NAN;
                    let x: $f = 0.5;
                    let y: $f = 1.0;

                    // Cases
                    let v = vec![neg_inf, x, neg_inf, neg_inf];
                    let rhs: $f = x;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v = vec![neg_inf, x, y, neg_inf];
                    let rhs: $f = y + ((x - y).exp() + 1.0).ln();
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v = vec![inf, x, y, neg_inf];
                    let rhs: $f = inf;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v = vec![x, inf, y, neg_inf];
                    let rhs: $f = inf;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v = vec![neg_inf, inf, x, y];
                    let rhs: $f = inf;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v = vec![neg_inf; 4];
                    let rhs: $f = neg_inf;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v = vec![inf; 4];
                    let rhs: $f = inf;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v = vec![neg_inf, neg_inf, neg_inf, inf];
                    let rhs: $f = inf;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v = vec![inf, neg_inf, neg_inf, neg_inf];
                    let rhs: $f = inf;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v = vec![inf, inf, y, neg_inf];
                    let rhs: $f = inf;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v = vec![inf, neg_inf, neg_inf, inf];
                    let rhs: $f = inf;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v = vec![neg_inf, neg_inf, neg_inf, x];
                    let rhs: $f = x;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v = vec![x, neg_inf, neg_inf, neg_inf];
                    let rhs: $f = x;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    // Edge cases
                    let v: Vec<$f> = vec![];
                    let rhs: $f = neg_inf;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v: Vec<$f> = vec![neg_inf];
                    let rhs: $f = neg_inf;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v: Vec<$f> = vec![inf];
                    let rhs: $f = inf;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v: Vec<$f> = vec![x];
                    let rhs: $f = x;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v: Vec<$f> = vec![0.0];
                    let rhs: $f = 0.0;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v: Vec<$f> = vec![x, inf];
                    let rhs: $f = inf;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    let v: Vec<$f> = vec![x, neg_inf];
                    let rhs: $f = x;
                    assert_eq!(v.iter().ln_sum_exp(), rhs);
                    assert_eq!(v.into_iter().ln_sum_exp(), rhs);

                    // Cases involving nan
                    let v = vec![x, inf, nan, y];
                    assert!(v.iter().ln_sum_exp().is_nan());
                    assert!(v.into_iter().ln_sum_exp().is_nan());

                    let v = vec![nan, x, y, inf];
                    assert!(v.iter().ln_sum_exp().is_nan());
                    assert!(v.into_iter().ln_sum_exp().is_nan());

                    let v = vec![inf, inf, x, nan];
                    assert!(v.iter().ln_sum_exp().is_nan());
                    assert!(v.into_iter().ln_sum_exp().is_nan());

                    let v = vec![inf, inf, neg_inf, nan];
                    assert!(v.iter().ln_sum_exp().is_nan());
                    assert!(v.into_iter().ln_sum_exp().is_nan());

                    let v = vec![inf, inf, inf, nan];
                    assert!(v.iter().ln_sum_exp().is_nan());
                    assert!(v.into_iter().ln_sum_exp().is_nan());

                    let v = vec![neg_inf, neg_inf, neg_inf, nan];
                    assert!(v.iter().ln_sum_exp().is_nan());
                    assert!(v.into_iter().ln_sum_exp().is_nan());

                    let v = vec![neg_inf, neg_inf, inf, nan];
                    assert!(v.iter().ln_sum_exp().is_nan());
                    assert!(v.into_iter().ln_sum_exp().is_nan());

                    let v = vec![x, y, nan, inf];
                    assert!(v.iter().ln_sum_exp().is_nan());
                    assert!(v.into_iter().ln_sum_exp().is_nan());

                    let v = vec![x, y, nan];
                    assert!(v.iter().ln_sum_exp().is_nan());
                    assert!(v.into_iter().ln_sum_exp().is_nan());

                    let v = vec![x, nan, y];
                    assert!(v.iter().ln_sum_exp().is_nan());
                    assert!(v.into_iter().ln_sum_exp().is_nan());

                    let v = vec![nan, x, y];
                    assert!(v.iter().ln_sum_exp().is_nan());
                    assert!(v.into_iter().ln_sum_exp().is_nan());

                    // Edge cases
                    let v: Vec<$f> = vec![nan];
                    assert!(v.iter().ln_sum_exp().is_nan());
                    assert!(v.into_iter().ln_sum_exp().is_nan());

                    let v: Vec<$f> = vec![nan; 4];
                    assert!(v.iter().ln_sum_exp().is_nan());
                    assert!(v.into_iter().ln_sum_exp().is_nan());

                    let v: Vec<$f> = vec![neg_inf, nan];
                    assert!(v.iter().ln_sum_exp().is_nan());
                    assert!(v.into_iter().ln_sum_exp().is_nan());

                    let v: Vec<$f> = vec![inf, nan];
                    assert!(v.iter().ln_sum_exp().is_nan());
                    assert!(v.into_iter().ln_sum_exp().is_nan());

                    let v: Vec<$f> = vec![neg_inf, inf, nan];
                    assert!(v.iter().ln_sum_exp().is_nan());
                    assert!(v.into_iter().ln_sum_exp().is_nan());

                    let v: Vec<$f> = vec![inf, neg_inf, nan];
                    assert!(v.iter().ln_sum_exp().is_nan());
                    assert!(v.into_iter().ln_sum_exp().is_nan());
                }

                #[test]
                fn ln_sum_exp_iterators_works() {
                    let v: Vec<$f> = vec![0.5, 1.0, 1.5];
                    let iter = v.iter().map(|&x| x);
                    assert_eq!(v.iter().ln_sum_exp(), iter.ln_sum_exp());

                    let iter = v.iter().map(|x| x.ln());
                    let rhs: $f = $f::ln(3.0);
                    assert!((iter.ln_sum_exp() -  rhs).abs() < 2.0 * $f::EPSILON);

                    let into_iter = v.into_iter().map(|x| x.ln());
                    let rhs: $f = $f::ln(3.0);
                    assert!((into_iter.ln_sum_exp() -  rhs).abs() < 2.0 * $f::EPSILON);

                    let into_iter = (1..4_i32).map(|x| x as $f).map(|x| x.ln());
                    let rhs: $f = $f::ln(6.0);
                    assert!((into_iter.ln_sum_exp() -  rhs).abs() < 2.0 * $f::EPSILON);

                    use std::collections::HashMap;
                    let map: HashMap<i32, $f> = HashMap::from([(1, 0.5), (2, 1.0), (3, 1.5)]);
                    let iter = map.values().map(|x| x.ln());
                    let rhs: $f = $f::ln(3.0);
                    assert!((iter.ln_sum_exp() -  rhs).abs() < 2.0 * $f::EPSILON);
                }
            }
        }
    }
    ln_sum_exp_tests! { f64_logsumexp_impl f64 }
    ln_sum_exp_tests! { f32_logsumexp_impl f32 }
}
