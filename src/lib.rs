use lnexp::LnExp;

pub trait LogAddExp<Rhs = Self> {
    type Output;
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

pub trait LogSumExp<T> {
    type Output;
    fn ln_sum_exp(self) -> Self::Output;
}

macro_rules! impl_logsumexp {
    { $($f:ident)+ } => {
        $(
            impl<U> LogSumExp<$f> for U
            where
                U: Iterator<Item = $f>,
            {
                type Output = $f;
                fn ln_sum_exp(mut self) -> Self::Output {
                    let mut m_old = $f::NEG_INFINITY;
                    let mut sum: $f = 0.0;
                    while let Some(v_i) = self.next() {
                        // if v_i != $f::NEG_INFINITY {
                        //     if v_i.is_nan() {
                        //         return v_i
                        //     } else {
                        //         let m_new = m_old.max(v_i);
                        //         sum = sum * (m_old - m_new).exp() + (v_i - m_new).exp();
                        //         m_old = m_new;
                        //     }
                        // }
                        if v_i != $f::NEG_INFINITY {
                            if v_i.is_nan() {
                                return v_i
                            } else if v_i == $f::INFINITY {
                                while let Some(v_i) = self.next() {
                                    if v_i.is_nan() {
                                        return v_i
                                    }
                                }
                                return $f::INFINITY
                            } else {
                                let m_new = m_old.max(v_i);
                                sum = sum * (m_old - m_new).exp() + (v_i - m_new).exp();
                                m_old = m_new;
                            }
                        }
                        // if v_i.is_finite() {
                        //     let m_new = m_old.max(v_i);
                        //     sum = sum * (m_old - m_new).exp() + (v_i - m_new).exp();
                        //     m_old = m_new;
                        // } else if v_i == $f::INFINITY {
                        //     while let Some(v_i) = self.next() {
                        //         if v_i.is_nan() {
                        //             return v_i
                        //         }
                        //     }
                        //     return $f::INFINITY
                        // } else if v_i.is_nan() {
                        //     return v_i
                        // }
                    }
                    m_old + sum.ln()
                    // if m_old == $f::INFINITY {
                    //     m_old
                    // } else {
                    //     m_old + sum.ln()
                    // }
                }
            }

            impl<'a, U> LogSumExp<&'a $f> for U
            where
                U: Iterator<Item = &'a $f>,
            {
                type Output = $f;
                fn ln_sum_exp(mut self) -> Self::Output {
                    let mut m_old = $f::NEG_INFINITY;
                    let mut sum: $f = 0.0;
                    while let Some(v_i) = self.next() {
                        if *v_i != $f::NEG_INFINITY {
                            if v_i.is_nan() {
                                return *v_i
                            } else {
                                let m_new = m_old.max(*v_i);
                                sum = sum * (m_old - m_new).exp() + (*v_i - m_new).exp();
                                m_old = m_new;
                            }
                        }
                    }
                    if m_old == $f::INFINITY {
                        m_old
                    } else {
                        m_old + sum.ln()
                    }
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

    #[cfg(test)]
    mod f64_logsumexp_impl {
        use super::*;

        #[test]
        fn ln_sum_exp_works() {
            let inf: f64 = f64::INFINITY;
            let neg_inf: f64 = f64::NEG_INFINITY;
            let nan: f64 = f64::NAN;
            let x: f64 = 0.5;
            let y: f64 = 1.0;

            // Cases
            let v = vec![neg_inf, x, neg_inf, neg_inf];
            let rhs: f64 = x;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v = vec![neg_inf, x, y, neg_inf];
            let rhs: f64 = y + ((x - y).exp() + 1.0).ln();
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v = vec![inf, x, y, neg_inf];
            let rhs: f64 = inf;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v = vec![x, inf, y, neg_inf];
            let rhs: f64 = inf;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v = vec![neg_inf, inf, x, y];
            let rhs: f64 = inf;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v = vec![neg_inf; 4];
            let rhs: f64 = neg_inf;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v = vec![inf; 4];
            let rhs: f64 = inf;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v = vec![neg_inf, neg_inf, neg_inf, inf];
            let rhs: f64 = inf;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v = vec![inf, neg_inf, neg_inf, neg_inf];
            let rhs: f64 = inf;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v = vec![inf, inf, y, neg_inf];
            let rhs: f64 = inf;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v = vec![inf, neg_inf, neg_inf, inf];
            let rhs: f64 = inf;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v = vec![neg_inf, neg_inf, neg_inf, x];
            let rhs: f64 = x;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v = vec![x, neg_inf, neg_inf, neg_inf];
            let rhs: f64 = x;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            // Edge cases
            let v: Vec<f64> = vec![];
            let rhs: f64 = neg_inf;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v: Vec<f64> = vec![neg_inf];
            let rhs: f64 = neg_inf;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v: Vec<f64> = vec![inf];
            let rhs: f64 = inf;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v: Vec<f64> = vec![x];
            let rhs: f64 = x;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v: Vec<f64> = vec![0.0];
            let rhs: f64 = 0.0;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v: Vec<f64> = vec![x, inf];
            let rhs: f64 = inf;
            assert_eq!(v.iter().ln_sum_exp(), rhs);
            assert_eq!(v.into_iter().ln_sum_exp(), rhs);

            let v: Vec<f64> = vec![x, neg_inf];
            let rhs: f64 = x;
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
            let v: Vec<f64> = vec![nan];
            assert!(v.iter().ln_sum_exp().is_nan());
            assert!(v.into_iter().ln_sum_exp().is_nan());

            let v: Vec<f64> = vec![nan; 4];
            assert!(v.iter().ln_sum_exp().is_nan());
            assert!(v.into_iter().ln_sum_exp().is_nan());

            let v: Vec<f64> = vec![neg_inf, nan];
            assert!(v.iter().ln_sum_exp().is_nan());
            assert!(v.into_iter().ln_sum_exp().is_nan());

            let v: Vec<f64> = vec![inf, nan];
            assert!(v.iter().ln_sum_exp().is_nan());
            assert!(v.into_iter().ln_sum_exp().is_nan());

            let v: Vec<f64> = vec![neg_inf, inf, nan];
            assert!(v.iter().ln_sum_exp().is_nan());
            assert!(v.into_iter().ln_sum_exp().is_nan());

            let v: Vec<f64> = vec![inf, neg_inf, nan];
            assert!(v.iter().ln_sum_exp().is_nan());
            assert!(v.into_iter().ln_sum_exp().is_nan());
        }
    }
}
