# Marginal and Conditional Multivariate Normal Distributions

Computes means and Cholesky factors of covariance or precision matrices
of multiple multivariate normal distributions.

## Usage

``` r
marg_mvnorm(chol, invchol, which = 1L)
cond_mvnorm(chol, invchol, which_given = 1L, given, center = FALSE)
```

## Arguments

- chol:

  Cholesky factors of covariance matrices as
  [`ltMatrices`](http://mvtnorm.R-forge.R-project.org/reference/ltMatrices.md)
  object, length is recylced to length of `lower` and `upper`.

- invchol:

  Cholesky factors of precision matrices as
  [`ltMatrices`](http://mvtnorm.R-forge.R-project.org/reference/ltMatrices.md)
  object, length is recylced to length of `lower` and `upper`. Either
  `chol` or `invchol` must be given.

- which:

  names or indices of elements those marginal distribution is of
  interest.

- which_given:

  names or indices of elements to condition on.

- given:

  matrix of realisations to condition on (number of rows is equal to
  `length(which)`, the number of columns corresponds to the number of
  matrices in `chol` or `invchol`.

- center:

  logical, if `TRUE`, the negative rescaled conditional mean is returned
  (such that it can be specified as `center` argument to
  [`slpmvnorm`](http://mvtnorm.R-forge.R-project.org/reference/lpmvnorm.md)).
  By default, the conditional mean is returned.

## Details

Derives parameters of the requested marginal or conditional
distributions, defined by `chol` (Cholesky factor of covariance) or
`invchol` (Cholesky factor of precision matrix) and, for conditional
distributions, the mean.

More details can be found in the `lmvnorm_src` package vignette.

## Value

A named list.

## See also

[`vignette("lmvnorm_src", package = "mvtnorm")`](https://cran.rstudio.com/web/packages/mvtnorm/vignettes/lmvnorm_src.pdf)
