# Multivariate Normal Density and Random Deviates

These functions provide the density function and a random number
generator for the multivariate normal distribution with mean equal to
`mean` and covariance matrix `sigma`.

## Usage

``` r
dmvnorm(x, mean = rep(0, p), sigma = diag(p), log = FALSE, checkSymmetry = TRUE)
rmvnorm(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
           method=c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, 
           checkSymmetry = TRUE, rnorm = stats::rnorm)
```

## Arguments

- x:

  vector or matrix of quantiles. When `x` is a matrix, each row is taken
  to be a quantile and columns correspond to the number of dimensions,
  `p`.

- n:

  number of observations.

- mean:

  mean vector, default is `rep(0, length = ncol(x))`. In `ldmvnorm` or
  `sldmvnorm`, `mean` is a matrix with observation-specific means
  arranged in columns.

- sigma:

  covariance matrix, default is `diag(ncol(x))`.

- log:

  logical; if `TRUE`, densities d are given as log(d).

- method:

  string specifying the matrix decomposition used to determine the
  matrix root of `sigma`. Possible methods are eigenvalue decomposition
  (`"eigen"`, default), singular value decomposition (`"svd"`), and
  Cholesky decomposition (`"chol"`). The Cholesky is typically fastest,
  not by much though.

- pre0.9_9994:

  logical; if `FALSE`, the output produced in mvtnorm versions up to
  0.9-9993 is reproduced. In 0.9-9994, the output is organized such that
  `rmvnorm(10,...)` has the same first ten rows as `rmvnorm(100, ...)`
  when called with the same seed.

- checkSymmetry:

  logical; if `FALSE`, skip checking whether the covariance matrix is
  symmetric or not. This will speed up the computation but may cause
  unexpected outputs when ill-behaved `sigma` is provided. The default
  value is `TRUE`.

- rnorm:

  a function with the same interface as
  [`rnorm`](https://rdrr.io/r/stats/Normal.html). This allows switching
  to other generators of standard normal variables.

## Details

`dmvnorm` computes the density function of the multivariate normal
specified by mean and the covariance matrix `sigma`.

`rmvnorm` generates multivariate normal variables.

## See also

[`pmvnorm`](http://mvtnorm.R-forge.R-project.org/reference/pmvnorm.md),
[`rnorm`](https://rdrr.io/r/stats/Normal.html),
[`qmvnorm`](http://mvtnorm.R-forge.R-project.org/reference/qmvnorm.md),
[`vignette("lmvnorm_src", package = "mvtnorm")`](https://cran.rstudio.com/web/packages/mvtnorm/vignettes/lmvnorm_src.pdf)

## Examples

``` r
dmvnorm(x=c(0,0))
#> [1] 0.1591549
dmvnorm(x=c(0,0), mean=c(1,1))
#> [1] 0.05854983

sigma <- matrix(c(4,2,2,3), ncol=2)
x <- rmvnorm(n=500, mean=c(1,2), sigma=sigma)
colMeans(x)
#> [1] 1.00517 2.01842
var(x)
#>          [,1]     [,2]
#> [1,] 4.068387 2.140247
#> [2,] 2.140247 3.085379
dS <- dmvnorm(x, sigma = sigma)

### alternative interface
C <- t(chol(sigma))
(C <- ltMatrices(C[lower.tri(C, diag = TRUE)], diag = TRUE))
#> , , 1
#> 
#>          1        2
#> 1 2.000000        .
#> 2 1.000000 1.414214
#> 
dC <- exp(ldmvnorm(obs = t(x), chol = C, logLik = FALSE))
all.equal(dS, dC)
#> [1] TRUE

x <- rmvnorm(n=500, mean=c(1,2), sigma=sigma, method="chol")
colMeans(x)
#> [1] 1.015820 1.951206
var(x)
#>          [,1]     [,2]
#> [1,] 3.676423 1.910980
#> [2,] 1.910980 2.898188

plot(x)
```
