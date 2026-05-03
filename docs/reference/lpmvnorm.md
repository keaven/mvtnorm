# Multivariate Normal Log-likelihood and Score Functions

Computes the log-likelihood (contributions) of multiple exact or
interval-censored observations (or a mix thereof) from multivariate
normal distributions and evaluates corresponding score functions.

## Usage

``` r
lpmvnorm(lower, upper, mean, invcholmean, center = NULL, chol, invchol, logLik = TRUE, 
         M = NULL, w = NULL, seed = NULL, tol = .Machine$double.eps, fast = FALSE)
slpmvnorm(lower, upper, mean, invcholmean, center = NULL, chol, invchol, logLik = TRUE, 
          M = NULL, w = NULL, seed = NULL, tol = .Machine$double.eps, fast = FALSE)
ldmvnorm(obs, mean, invcholmean, chol, invchol, logLik = TRUE) 
sldmvnorm(obs, mean, invcholmean, chol, invchol, logLik = TRUE) 
ldpmvnorm(obs, lower, upper, mean, invcholmean, chol, invchol, logLik = TRUE, ...) 
sldpmvnorm(obs, lower, upper, mean, invcholmean, chol, invchol, logLik = TRUE, ...)
```

## Arguments

- lower:

  matrix of lower limits (one column for each observation, \\J\\ rows).

- upper:

  matrix of upper limits (one column for each observation, \\J\\ rows).

- obs:

  matrix of exact observations (one column for each observation, \\J\\
  rows).

- mean:

  matrix of means (one column for each observation, length is recycled
  to length of `obs`, `lower` and `upper`).

- invcholmean:

  matrix of means left-multiplied with inverse Cholesky factor
  (`invchol %*% mean`, one column for each observation, length is
  recycled to length of `obs`, `lower` and `upper`).

- center:

  matrix of negative rescaled means (one column for each observation,
  length is recycled to length of `lower` and `upper`) as returned by
  `cond_mvnorm(..., center = TRUE).`.

- chol:

  Cholesky factors of covariance matrices as
  [`ltMatrices`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  object, length is recylced to length of `obs`, `lower` and `upper`.

- invchol:

  Cholesky factors of precision matrices as
  [`ltMatrices`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  object, length is recylced to length of `lower` and `upper`. Either
  `chol` or `invchol` must be given.

- logLik:

  logical, if `TRUE`, the log-likelihood is returned, otherwise the
  individual contributions to the sum are returned.

- M:

  number of iterations, early stopping based on estimated errors is NOT
  implemented.

- w:

  an optional matrix of weights with \\J - 1\\ rows. This allows to
  replace the default Monte-Carlo procedure mvtnorm::numerical-:1992
  with a quasi-Monte-Carlo approach mvtnorm::Genz_Bretz_2002. Note that
  the same weights for evaluating the multivariate normal probability
  are used for all observations when `ncol(w) == M` is specified. If
  `ncol(w) == ncol(lower) * M`, each likelihood contribution is
  evaluated on the corresponding sub-matrix. If `w` is `NULL`, different
  uniform numbers are drawn for each observation.

- seed:

  an object specifying if and how the random number generator should be
  initialized, see [`simulate`](https://rdrr.io/r/stats/simulate.html).
  Only applied when `w` is `NULL`.

- tol:

  tolerance limit, values smaller than `tol` are interpreted as zero.

- fast:

  logical, if `TRUE`, a faster but less accurate version of `pnorm` is
  used internally.

- ...:

  additional arguments to `lpmvnorm`.

## Details

Evaluates the multivariate normal log-likelihood defined by `means` and
`chol` over boxes defined by `lower` and `upper` or for exact
observations `obs`.

Monte-Carlo \|mvtnorm::numerical-:1992\|the default and
quasi-Monte-Carlo mvtnorm::Genz_Bretz_2002 integration is implemented,
the latter with weights obtained, for example, from packages qrng or
randtoolbox. It is the responsibility of the user to ensure a meaningful
lattice is used. In case of doubt, use plain Monte-Carlo (`w = NULL`) or
[`pmvnorm`](https://keaven.github.io/mvtnorm/reference/pmvnorm.md).

`slpmvnorm` computes both the individual log-likelihood contributions
and the corresponding score matrix (of dimension \\J \times (J + 1) / 2
\times N\\) if `chol` contains diagonal elements. Otherwise, the
dimension is \\J \times (J - 1) / 2 \times N\\. The scores for exact or
mixed exact-interval observations are computed by `sldmvnorm` and
`sldpmvnorm`, respectively.

More details can be found in the `lmvnorm_src` package vignette.

## Value

The log-likelihood (`logLik = TRUE`) or the individual contributions to
the log-likelihood. `slpmvnorm`, `sldmvnorm`, and `sldpmvnorm` return
the score matrices and, optionally (`logLik = TRUE`), the individual
log-likelihood contributions as well as scores for `obs`, `lower`,
`upper`, and `mean`.

## See also

[`dmvnorm`](https://keaven.github.io/mvtnorm/reference/Mvnorm.md),
[`vignette("lmvnorm_src", package = "mvtnorm")`](https://cran.rstudio.com/web/packages/mvtnorm/vignettes/lmvnorm_src.pdf)

## References

\*

## Examples

``` r
  ### five observations
  N <- 5L
  ### dimension
  J <- 4L

  ### lower and upper bounds, ie interval-censoring
  lwr <- matrix(-runif(N * J), nrow = J)
  upr <- matrix(runif(N * J), nrow = J)

  ### Cholesky factor
  (C <- ltMatrices(runif(J * (J + 1) / 2), diag = TRUE))
#> , , 1
#> 
#>            1          2          3          4
#> 1 0.04904382          .          .          .
#> 2 0.36110032 0.96296355          .          .
#> 3 0.99842188 0.08205368 0.44223926          .
#> 4 0.10186426 0.58576477 0.66234805 0.41829153
#> 
  ### corresponding covariance matrix
  (S <- as.array(Tcrossprod(C))[,,1])
#>             1          2          3           4
#> 1 0.002405297 0.01770974 0.04896643 0.004995813
#> 2 0.017709740 1.05769224 0.43954516 0.600853333
#> 3 0.048966425 0.43954516 1.19915463 0.442683968
#> 4 0.004995813 0.60085333 0.44268397 0.967169423

  ### plain Monte-Carlo (Genz, 1992)
  w <- NULL
  M <- 25000
  ### quasi-Monte-Carlo (Genz & Bretz, 2002, but with different weights)
  if (require("qrng")) w <- t(ghalton(M * N, J - 1))
#> Loading required package: qrng

  ### log-likelihood
  lpmvnorm(lower = lwr, upper = upr, chol = C, w = w, M = M)
#> [1] -16.78885

  ### compare with pmvnorm
  exp(lpmvnorm(lower = lwr, upper = upr, chol = C, logLik = FALSE, w = w, M = M))
#> [1] 0.02322384 0.04121074 0.04257901 0.03160109 0.03970588
  sapply(1:N, function(i) pmvnorm(lower = lwr[,i], upper = upr[,i], sigma = S))
#> [1] 0.02321819 0.04121074 0.04257357 0.03160343 0.03970454

  ### log-lik contributions and score matrix
  slpmvnorm(lower = lwr, upper = upr, chol = C, w = w, M = M, logLik = TRUE)
#> $logLik
#> [1] -3.762576 -3.189056 -3.156394 -3.454564 -3.226256
#> 
#> $mean
#>         [,1]          [,2]         [,3]          [,4]          [,5]
#> 1 -2.5632952  1.111698e-10  0.001107593  3.471900e-32  8.637699e-05
#> 2  0.4311371  3.263094e-01 -0.218912451  4.540177e-01 -1.987211e-01
#> 3  0.4658949 -8.394498e-02 -0.010388578  3.490625e-01  3.691141e-01
#> 4 -0.5561359 -1.924189e-01  0.285118334 -7.239072e-01 -2.177145e-01
#> 
#> $lower
#>            [,1]          [,2]         [,3]          [,4]          [,5]
#> 1  8.493698e-15 -1.111698e-10 -0.001107593 -8.131617e-62 -9.890908e-05
#> 2 -9.109586e-01 -8.845568e-01 -0.449826839 -1.686904e+00 -1.624378e+00
#> 3 -1.072127e+00 -2.272837e+00 -0.756021043 -1.331178e+00 -1.400678e+00
#> 4 -2.866038e+00 -5.093898e-01 -2.560435964 -5.029333e-01 -3.769746e-01
#> 
#> $upper
#>        [,1]          [,2]         [,3]          [,4]         [,5]
#> 1 2.5632952 -1.650267e-80 4.719704e-57 -3.471900e-32 0.0000125321
#> 2 0.4798215  5.582474e-01 6.687393e-01  1.232887e+00 1.8230992263
#> 3 0.6062318  2.356782e+00 7.664096e-01  9.821158e-01 1.0315640081
#> 4 3.4221741  7.018087e-01 2.275318e+00  1.226841e+00 0.5946891633
#> 
#> $chol
#> , , 1
#> 
#>               1             2             3             4
#> 1 -3.128377e+00             .             .             .
#> 2  1.538494e-02 -7.357650e-01             .             .
#> 3 -5.503310e-01  3.350605e-01 -3.718220e-01             .
#> 4  1.113629e-01 -3.107618e-01 -7.838638e-01 -5.936729e-01
#> 
#> , , 2
#> 
#>               1             2             3             4
#> 1 -6.507157e-10             .             .             .
#> 2 -1.385765e-01 -7.591018e-01             .             .
#> 3 -8.935289e-01  2.863887e-01 -2.544427e-01             .
#> 4  3.321968e-01 -2.366361e-01 -6.969062e-01 -5.195500e-01
#> 
#> , , 3
#> 
#>               1             2             3             4
#> 1 -1.808960e-03             .             .             .
#> 2 -1.687015e-01 -7.291574e-01             .             .
#> 3 -8.212215e-01  2.735014e-01 -2.048929e-01             .
#> 4  3.611778e-01 -2.804092e-01 -8.369607e-01 -6.142420e-01
#> 
#> , , 4
#> 
#>               1             2             3             4
#> 1  4.045095e-31             .             .             .
#> 2  1.309706e-02 -8.967016e-01             .             .
#> 3 -7.045883e-01  3.270607e-01 -3.676760e-01             .
#> 4  4.806793e-02 -1.935461e-01 -5.960965e-01 -3.800637e-01
#> 
#> , , 5
#> 
#>               1             2             3             4
#> 1 -1.571018e-04             .             .             .
#> 2 -2.144170e-01 -8.836940e-01             .             .
#> 3 -7.495461e-01  2.050522e-01 -2.746493e-01             .
#> 4  2.501904e-01 -6.944593e-02 -7.062036e-01 -4.863555e-01
#> 
#> 
```
