# Quantiles of the Multivariate t Distribution

Computes the equicoordinate quantile function of the multivariate t
distribution for arbitrary correlation matrices based on inversion of
[`pmvt`](https://keaven.github.io/mvtnorm/reference/pmvt.md), using a
stochastic root finding algorithm described in Bornkamp (2018).

## Usage

``` r
qmvt(p, interval = NULL, tail = c("lower.tail", "upper.tail", "both.tails"), 
     df = 1, delta = 0, corr = NULL, sigma = NULL, algorithm = GenzBretz(),
     type = c("Kshirsagar", "shifted"), ptol = 0.001, maxiter = 500, 
     trace = FALSE, seed = NULL, ...)
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

- delta:

  the vector of noncentrality parameters of length n, for
  `type = "shifted"` delta specifies the mode.

- df:

  degree of freedom as integer. Normal quantiles are computed for
  `df = 0` or `df = Inf`.

- corr:

  the correlation matrix of dimension n.

- sigma:

  the covariance matrix of dimension n. Either `corr` or `sigma` can be
  specified. If `sigma` is given, the problem is standardized
  internally. If `corr` is given, it is assumed that appropriate
  standardization was performed by the user. If neither `corr` nor
  `sigma` is given, the identity matrix in the univariate case (so
  `corr = 1`) is used for `corr`.

- algorithm:

  an object of class
  [`GenzBretz`](https://keaven.github.io/mvtnorm/reference/algorithms.md)
  or
  [`TVPACK`](https://keaven.github.io/mvtnorm/reference/algorithms.md)
  defining the hyper parameters of this algorithm.

- type:

  type of the noncentral multivariate t distribution to be computed. The
  choice `type = "Kshirsagar"` corresponds to formula (1.4) in
  mvtnorm::Genz_Bretz_2009, see also Chapter 5.1 in
  mvtnorm::Kotz+Nadarajah:2004. This is the noncentral t-distribution
  needed for calculating the power of multiple contrast tests under a
  normality assumption. `type = "shifted"` corresponds to the formula
  right before formula (1.4) in mvtnorm::Genz_Bretz_2009 see also
  formula (1.1) in\|mvtnorm::Kotz+Nadarajah:2004\|.

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
[`qmvnorm`](https://keaven.github.io/mvtnorm/reference/qmvnorm.md)

## Examples

``` r
## basic evaluation
qmvt(0.95, df = 16, tail = "both")
#> $quantile
#> [1] 2.119905
#> 
#> $f.quantile
#> [1] 0.975
#> 

## check behavior for df=0 and df=Inf
Sigma <- diag(2)
set.seed(29)
q0 <- qmvt(0.95, sigma = Sigma, df = 0,   tail = "both")$quantile
set.seed(29)
q8 <- qmvt(0.95, sigma = Sigma, df = Inf, tail = "both")$quantile
set.seed(29)
qn <- qmvnorm(0.95, sigma = Sigma, tail = "both")$quantile
stopifnot(identical(q0, q8),
          isTRUE(all.equal(q0, qn, tol = (.Machine$double.eps)^(1/3))))

## if neither sigma nor corr are provided, corr = 1 is used internally
df <- 0
set.seed(29)
qt95 <- qmvt(0.95, df = df, tail = "both")$quantile
set.seed(29)
qt95.c <- qmvt(0.95, df = df, corr  = 1, tail = "both")$quantile
set.seed(29)
qt95.s <- qmvt(0.95, df = df, sigma = 1, tail = "both")$quantile
stopifnot(identical(qt95, qt95.c),
          identical(qt95, qt95.s))

df <- 4
set.seed(29)
qt95 <- qmvt(0.95, df = df, tail = "both")$quantile
set.seed(29)
qt95.c <- qmvt(0.95, df = df, corr  = 1, tail = "both")$quantile
set.seed(29)
qt95.s <- qmvt(0.95, df = df, sigma = 1, tail = "both")$quantile
stopifnot(identical(qt95, qt95.c),
          identical(qt95, qt95.s))
```
