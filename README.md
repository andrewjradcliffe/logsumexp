# logsumexp

## Description

Provides two traits:
* `LogAddExp`, which provides a numerically stable evaluation of
  `log(exp(a) + exp(b))`, utilizing `ln_1p_exp` from the `lnexp`
  crate. Implementations are provided for for `f64` and `f32` (and
  their respective references).
* `LogSumExp`, which provides a numerically stable, 1-pass algorithm
  for evaluation of
  [LogSumExp](https://en.wikipedia.org/wiki/LogSumExp) with correct
  handling of +/- infinity and `nan`. Implementations are provided
  iterators which produce `Item`s of `f64` or `f32` (and for
  respective references).

## License

Licensed under either of

  * [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0)
  * [MIT license](http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.

## Citations
- [Milakov, Maxim, and Natalia Gimelshein. "Online normalizer calculation for softmax." arXiv preprint arXiv:1805.02867 (2018)](https://arxiv.org/pdf/1805.02867.pdf)
