use lnexp::LnExp;
pub trait LogAddExp {
    fn ln_add_exp(&self, other: Self) -> Self;
}

macro_rules! impl_logaddexp {
    ( $f:ident ) => {
        impl LogAddExp for $f {
            fn ln_add_exp(&self, other: $f) -> $f {
                // if *self == other {
                //     *self + std::$f::consts::LN_2
                // } else {
                //     let (max, diff) = if *self < other {
                //         (other, *self - other)
                //     } else {
                //         (*self, other - *self)
                //     };
                //     max + diff.ln_1p_exp()
                // }
                let (max, diff) = if *self < other {
                    (other, *self - other)
                } else {
                    if *self == other {
                        (other, 0.0)
                    } else {
                        (*self, other - *self)
                    }
                };
                max + diff.ln_1p_exp()
            }
        }
        // impl<'a> LogAddExp for &'a $f {
        //     fn ln_add_exp(&self, other: &'a $f) -> $f {
        //         self.ln_add_exp(*other)
        //     }
        // }
    };
}
impl_logaddexp!(f64);
impl_logaddexp!(f32);

// pub trait LogSumExp<A = Self>: Sized {
//     fn ln_sum_exp<I: Iterator<Item = A>>(iter: I) -> Self;
// }

// impl LogSumExp for f64 {
//     fn ln_sum_exp<I: Iterator<Item = Self>>(iter: I) -> Self {
//         let (m_old, sum) = iter.fold((f64::NEG_INFINITY, 0.0), |(m_old, sum), v_i| {
//             if v_i != f64::NEG_INFINITY {
//                 let m_new = m_old.max(v_i);
//                 (m_new, sum * (m_old - m_new).exp() + (v_i - m_new).exp())
//             } else {
//                 (m_old, sum)
//             }
//         });
//         if m_old == f64::INFINITY {
//             m_old
//         } else {
//             m_old + sum.ln()
//         }
//     }
// }

// pub trait LogSumExp {
//     type Output;

//     fn ln_sum_exp(self) -> Self::Output;
// }

// impl<T> LogSumExp for T
// where
//     T: Iterator<Item = f64>,
// {
//     type Output = f64;
//     fn ln_sum_exp(self) -> Self::Output {
//         let (m_old, sum) = self.fold((f64::NEG_INFINITY, 0.0), |(m_old, sum), v_i| {
//             if v_i != f64::NEG_INFINITY {
//                 let m_new = m_old.max(v_i);
//                 (m_new, sum * (m_old - m_new).exp() + (v_i - m_new).exp())
//             } else {
//                 (m_old, sum)
//             }
//         });
//         if m_old == f64::INFINITY {
//             m_old
//         } else {
//             m_old + sum.ln()
//         }
//     }
// }

// impl<'a, T> LogSumExp for T
// where
//     T: Iterator<Item = &'a f64>,
// {
//     type Output = f64;
//     fn ln_sum_exp(self) -> Self::Output {
//         let (m_old, sum) = self.fold((f64::NEG_INFINITY, 0.0), |(m_old, sum), v_i| {
//             if *v_i != f64::NEG_INFINITY {
//                 let m_new = m_old.max(*v_i);
//                 (m_new, sum * (m_old - m_new).exp() + (*v_i - m_new).exp())
//             } else {
//                 (m_old, sum)
//             }
//         });
//         if m_old == f64::INFINITY {
//             m_old
//         } else {
//             m_old + sum.ln()
//         }
//     }
// }

pub trait LogSumExp<T> {
    fn ln_sum_exp(self) -> T;
}
impl<U> LogSumExp<f64> for U
where
    U: Iterator<Item = f64>,
{
    // type Output = f64;
    fn ln_sum_exp(self) -> f64 {
        let (m_old, sum) = self.fold((f64::NEG_INFINITY, 0.0), |(m_old, sum), v_i| {
            if v_i != f64::NEG_INFINITY {
                let m_new = m_old.max(v_i);
                (m_new, sum * (m_old - m_new).exp() + (v_i - m_new).exp())
            } else {
                (m_old, sum)
            }
        });
        if m_old == f64::INFINITY {
            m_old
        } else {
            m_old + sum.ln()
        }
    }
}

// impl<'a, T> LogSumExp<'a> for T
// where
//     T: Iterator<Item = &'a f64>,
// {
//     type Output = f64;
//     fn ln_sum_exp(self) -> Self::Output {
//         let (m_old, sum) = self.fold((f64::NEG_INFINITY, 0.0), |(m_old, sum), v_i| {
//             if *v_i != f64::NEG_INFINITY {
//                 let m_new = m_old.max(*v_i);
//                 (m_new, sum * (m_old - m_new).exp() + (*v_i - m_new).exp())
//             } else {
//                 (m_old, sum)
//             }
//         });
//         if m_old == f64::INFINITY {
//             m_old
//         } else {
//             m_old + sum.ln()
//         }
//     }
// }

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
