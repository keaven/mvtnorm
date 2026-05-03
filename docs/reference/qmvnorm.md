# Quantiles of the Multivariate Normal Distribution

Computes the equicoordinate quantile function of the multivariate normal
distribution for arbitrary correlation matrices based on inversion of
[`pmvnorm`](https://keaven.github.io/mvtnorm/reference/pmvnorm.md),
using a stochastic root finding algorithm described in Bornkamp (2018).

## Usage

``` r
qmvnorm(p, interval = NULL, tail = c("lower.tail", "upper.tail", "both.tails"), 
        mean = 0, corr = NULL, sigma = NULL, algorithm = GenzBretz(),
        ptol = 0.001, maxiter = 500, trace = FALSE, seed = NULL, ...)
```

## Arguments

- p:

  probability.

- interval:

  optional, a vector containing the end-points of the interval to be
  searched. Does not need to contain the true quantile, just used as
  starting values by the root-finder. If equal to NULL a guess is used.

- tail:

  specifies which quantiles should be computed. `lower.tail` gives the
  quantile \\x\\ for which \\P\[X \le x\] = p\\, `upper.tail` gives
  \\x\\ with \\P\[X \> x\] = p\\ and `both.tails` leads to \\x\\ with
  \\P\[-x \le X \le x\] = p\\.

- mean:

  the mean vector of length n.

- corr:

  the correlation matrix of dimension n.

- sigma:

  the covariance matrix of dimension n. Either `corr` or `sigma` can be
  specified. If `sigma` is given, the problem is standardized
  internally. If `corr` is given, it is assumed that appropriate
  standardization was performed by the user. If neither `corr` nor
  `sigma` is given, the identity matrix is used for `sigma`.

- algorithm:

  an object of class
  [`GenzBretz`](https://keaven.github.io/mvtnorm/reference/algorithms.md),
  [`Miwa`](https://keaven.github.io/mvtnorm/reference/algorithms.md) or
  [`TVPACK`](https://keaven.github.io/mvtnorm/reference/algorithms.md)
  specifying both the algorithm to be used as well as the associated
  hyper parameters.

- ptol, maxiter, trace:

  Parameters passed to the stochastic root-finding algorithm. Iteration
  stops when the 95% confidence interval for the predicted quantile is
  inside \[p-ptol, p+ptol\]. `maxiter` is the maximum number of
  iterations for the root finding algorithm. `trace` prints the
  iterations of the root finder.

- seed:

  an object specifying if and how the random number generator should be
  initialized, see [`simulate`](https://rdrr.io/r/stats/simulate.html).

- ...:

  additional parameters to be passed to
  [`GenzBretz`](https://keaven.github.io/mvtnorm/reference/algorithms.md).

## Details

Only equicoordinate quantiles are computed, i.e., the quantiles in each
dimension coincide. The result is seed dependend. The concept is
explained in mvtnorm::Bornkamp:2018.

## Value

A list with two components: `quantile` and `f.quantile` give the
location of the quantile and the difference between the distribution
function evaluated at the quantile and `p`.

## References

\*

## See also

[`pmvnorm`](https://keaven.github.io/mvtnorm/reference/pmvnorm.md),
[`qmvt`](https://keaven.github.io/mvtnorm/reference/qmvt.md)

## Examples

``` r
qmvnorm(0.95, sigma = diag(2), tail = "both")
#> $quantile
#> [1] 2.236422
#> 
#> $f.quantile
#> [1] -1.310417e-06
#> 
#> attr(,"message")
#> [1] "Normal Completion"
```
