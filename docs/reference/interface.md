# User Interface to Multiple Multivariate Normal Distributions

A simple user interface for computing on multiple multivariate normal
distributions.

## Usage

``` r
mvnorm(mean, invcholmean, chol, invchol)
# S3 method for class 'mvnorm'
aperm(a, perm, ...)
margDist(object, which, ...)
# S3 method for class 'mvnorm'
margDist(object, which, ...)
condDist(object, which_given, given, ...)
# S3 method for class 'mvnorm'
condDist(object, which_given, given, ...)
# S3 method for class 'mvnorm'
simulate(object, nsim = dim(object$scale)[1L], seed = NULL, 
                            standardize = FALSE, as.data.frame = FALSE, ...)
# S3 method for class 'mvnorm'
logLik(object, obs, lower, upper, standardize = FALSE, ...)
# S3 method for class 'mvnorm'
lLgrad(object, obs, lower, upper, standardize = FALSE, ...)
```

## Arguments

- chol:

  either an `ltMatrices` object specifying (multiple) Cholesky factors
  of the covariance matrix or one single numeric lower triangular square
  matrix.

- invchol:

  either an `ltMatrices` object specifying (multiple) inverse Cholesky
  factors of the covariance matrix or one single numeric lower
  triangular square matrix.

- a,object:

  objects of class `mvnorm`.

- perm:

  a permutation of the covariance matrix corresponding to `a`.

- which:

  names or indices of elements those marginal distribution is of
  interest.

- which_given:

  names or indices of elements to condition on.

- given:

  matrix of realisations to condition on (number of rows is equal to
  `length(which)`, the number of columns corresponds to the number of
  matrices in `chol` or `invchol`.

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

  matrix of scaled means, that is, `invchol %*% mean`, one column for
  each observation, length is recycled to length of `obs`, `lower` and
  `upper`.

- seed:

  an object specifying if and how the random number generator should be
  initialized, see [`simulate`](https://rdrr.io/r/stats/simulate.html).

- standardize:

  logical, should the Cholesky factor (or its inverse) undergo
  standardization (ensuring the covariance matrix is a correlation
  matrix) before computing the likelihood.

- nsim:

  number of samples to draw.

- as.data.frame:

  logical, convert the \$J x N\$ matrix result to a classical \$N x J\$
  data frame.

- ...:

  Additional arguments to
  [`ldpmvnorm`](https://keaven.github.io/mvtnorm/reference/lpmvnorm.md)
  and
  [`sldpmvnorm`](https://keaven.github.io/mvtnorm/reference/lpmvnorm.md)

## Details

The constructor `mvnorm` can be used to specify (multiple) multivariate
normal distributions. `margDist` derives marginal and `condDist`
conditional distributions from such objects. A `simulate` method exists
for drawn samples from multivariate normals.

The continuous (data in `obs`), discrete (intervals in `lower` and
`upper`), and mixed continuous-discrete log-likelihood is implemented in
`logLik`. The corresponding gradients with respect to all model
parameters and with respect to the data arguments is available from
`lLgrad`.

Rationals and examples are given in Chapter 7 of the package vignette
linked to below.

## Value

`mvnorm`, `margDist`, and `condDist` return objects of class `mvnorm`.
`logLik` returns the log-likelihood and `lLgrad` a list with gradients.

## See also

[`vignette("lmvnorm_src", package = "mvtnorm")`](https://cran.rstudio.com/web/packages/mvtnorm/vignettes/lmvnorm_src.pdf)
