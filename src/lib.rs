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
                        if v_i != $f::NEG_INFINITY {
                            if v_i.is_nan() {
                                return v_i
                            } else {
                                let m_new = m_old.max(v_i);
                                sum = sum * (m_old - m_new).exp() + (v_i - m_new).exp();
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
                    let nan: $f = 0.0 / 0.0;
                    let x: $f = inf;
                    let y: $f = 0.5;
                    assert_eq!(x.ln_add_exp(inf), inf);
                    assert_eq!(x.ln_add_exp(-inf), inf);
                    assert_eq!(x.ln_add_exp(y), inf);
                    assert_eq!(y.ln_add_exp(x), inf);

                    let x: $f = -inf;
                    assert_eq!(x.ln_add_exp(inf), inf);
                    assert_eq!(x.ln_add_exp(-inf), -inf);
                    assert_eq!(x.ln_add_exp(y), y);
                    assert_eq!(y.ln_add_exp(x), y);

                    assert!(nan.ln_add_exp(inf).is_nan());
                    assert!(nan.ln_add_exp(-inf).is_nan());
                    assert!(inf.ln_add_exp(nan).is_nan());
                    assert!((-inf).ln_add_exp(nan).is_nan());

                    let x: $f = 0.65;
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
            assert_eq!(v1.into_iter().ln_sum_exp(), f64::NEG_INFINITY);

            let v2 = vec![f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY, 0.5];
            assert_eq!(v2.into_iter().ln_sum_exp(), 0.5);

            let v3 = vec![0.5, f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY];
            assert_eq!(v3.into_iter().ln_sum_exp(), 0.5);

            let v4 = vec![
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
            ];
            assert_eq!(v4.into_iter().ln_sum_exp(), f64::INFINITY);

            let v5 = vec![f64::INFINITY; 4];
            assert_eq!(v5.into_iter().ln_sum_exp(), f64::INFINITY);

            let v6 = vec![
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
            ];
            assert_eq!(v6.into_iter().ln_sum_exp(), f64::INFINITY);
        }
    }
}
