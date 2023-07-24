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
            let v1 = vec![f64::NEG_INFINITY; 4];
            assert_eq!(v1.iter().ln_sum_exp(), f64::NEG_INFINITY);
            assert_eq!(v1.into_iter().ln_sum_exp(), f64::NEG_INFINITY);

            let v2 = vec![f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY, 0.5];
            assert_eq!(v2.iter().ln_sum_exp(), 0.5);
            assert_eq!(v2.into_iter().ln_sum_exp(), 0.5);

            let v3 = vec![0.5, f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY];
            assert_eq!(v3.iter().ln_sum_exp(), 0.5);
            assert_eq!(v3.into_iter().ln_sum_exp(), 0.5);

            let v4 = vec![
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
            ];
            assert_eq!(v4.iter().ln_sum_exp(), f64::INFINITY);
            assert_eq!(v4.into_iter().ln_sum_exp(), f64::INFINITY);

            let v5 = vec![f64::INFINITY; 4];
            assert_eq!(v5.iter().ln_sum_exp(), f64::INFINITY);
            assert_eq!(v5.into_iter().ln_sum_exp(), f64::INFINITY);

            let v6 = vec![
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
            ];
            assert_eq!(v6.iter().ln_sum_exp(), f64::INFINITY);
            assert_eq!(v6.into_iter().ln_sum_exp(), f64::INFINITY);

            let v7: Vec<f64> = vec![];
            assert_eq!(v7.iter().ln_sum_exp(), f64::NEG_INFINITY);
            assert_eq!(v7.into_iter().ln_sum_exp(), f64::NEG_INFINITY);

            let v8 = vec![
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
            ];
            assert_eq!(v8.iter().ln_sum_exp(), f64::INFINITY);
            assert_eq!(v8.into_iter().ln_sum_exp(), f64::INFINITY);

            let nan: f64 = 0.0 / 0.0;
            let v9 = vec![1.5, 1.0, nan];
            assert!(v9.iter().ln_sum_exp().is_nan());
            assert!(v9.into_iter().ln_sum_exp().is_nan());
            let v10 = vec![nan];
            assert!(v10.iter().ln_sum_exp().is_nan());
            assert!(v10.into_iter().ln_sum_exp().is_nan());
            let v11 = vec![nan, 1.5, 1.0];
            assert!(v11.iter().ln_sum_exp().is_nan());
            assert!(v11.into_iter().ln_sum_exp().is_nan());
            let v12 = vec![1.5, nan, 1.0];
            assert!(v12.iter().ln_sum_exp().is_nan());
            assert!(v12.into_iter().ln_sum_exp().is_nan());

            let v13 = vec![nan; 3];
            assert!(v13.iter().ln_sum_exp().is_nan());
            assert!(v13.into_iter().ln_sum_exp().is_nan());

            let v14 = vec![nan, 1.5, f64::INFINITY];
            assert!(v14.iter().ln_sum_exp().is_nan());
            assert!(v14.into_iter().ln_sum_exp().is_nan());

            let v15 = vec![nan, f64::INFINITY];
            assert!(v15.iter().ln_sum_exp().is_nan());
            assert!(v15.into_iter().ln_sum_exp().is_nan());

            let v16 = vec![f64::INFINITY, f64::NEG_INFINITY, nan];
            assert!(v16.iter().ln_sum_exp().is_nan());
            assert!(v16.into_iter().ln_sum_exp().is_nan());

            let v17 = vec![f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, nan];
            assert!(v17.iter().ln_sum_exp().is_nan());
            assert!(v17.into_iter().ln_sum_exp().is_nan());
        }
    }
}
