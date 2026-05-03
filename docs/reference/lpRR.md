# Multivariate Normal Log-likelihood and Score Functions for Reduced Rank Covariances

Computes the log-likelihood (contributions) of interval-censored
observations from multivariate normal distributions with reduced rank
structure and evaluates corresponding score functions.

## Usage

``` r
lpRR(lower, upper, mean = 0, B, D = rep(1, nrow(B)), 
     Z, weights = 1 / ncol(Z), log.p = TRUE)
slpRR(lower, upper, mean = 0, B, D = rep(1, nrow(B)), 
      Z, weights = 1 / ncol(Z), log.p = TRUE)
```

## Arguments

- lower:

  vector of lower limits (one element for each dimension, \\J\\
  elements).

- upper:

  vector of upper limits (one element for each dimension, \\J\\
  elements).

- mean:

  vector of means (one element for each dimension, length is recycled to
  length of `lower` and `upper`).

- B:

  matrix of dimension \\J \times K\\.

- D:

  vector of \\J\\ diagonal elements.

- Z:

  matrix of standard normal random variables, with \\K\\ nrows.

- weights:

  optional weights.

- log.p:

  logical. By default, log-probabilities are returned.

## Details

Evaluates the multivariate normal log-likelihood defined by `mean`, `B`
and `D` when the covariance is \\\Sigma = B B^\top + D\\ over boxes
defined by `lower` and `upper`. Details are given in
mvtnorm::Genz_Bretz_2009, Chapter 2.3.1.

`slpRR` computes the corresponding score functions with respect to
`lower`, `upper`, `mean`, `B` and `D`.

More details can be found in the `lmvnorm_src` package vignette.

## Value

The log-likelihood (`log.p = TRUE`) or corresponding probability.
`slpRR` return the scores.

## See also

[`vignette("lmvnorm_src", package = "mvtnorm")`](https://cran.rstudio.com/web/packages/mvtnorm/vignettes/lmvnorm_src.pdf)

## References

\*

## Examples

``` r
  J <- 6
  K <- 3
  B <- matrix(rnorm(J * K), nrow = J)
  D <- runif(J)
  S <- tcrossprod(B) + diag(D)
  a <- -(2 + runif(J))
  b <- 2 + runif(J)
  M <- 1e4
  Z <- matrix(rnorm(K * M), nrow = K)
  ## log-likelihood
  lpRR(lower = a, upper = b, B = B, D = D, Z = Z)
#> [1] -0.8808033
  ## score wrt all arguments
  slpRR(lower = a, upper = b, B = B, D = D, Z = Z)
#> $lower
#> [1] -0.172229168 -0.095932619 -0.005381103 -0.042516356 -0.031346735
#> [6] -0.060970576
#> 
#> $upper
#> [1] 0.169332988 0.109853513 0.005504726 0.010142701 0.011226195 0.057602701
#> 
#> $mean
#> [1]  0.0028961792 -0.0139208948 -0.0001236227  0.0323736549  0.0201205403
#> [6]  0.0033678752
#> 
#> $B
#>              [,1]         [,2]          [,3]
#> [1,] -0.190824059  0.181469511 -0.0218882493
#> [2,]  0.085340496  0.058709862 -0.1835994404
#> [3,] -0.005001726 -0.010939550 -0.0168299430
#> [4,]  0.047841993  0.026602709 -0.0426801439
#> [5,] -0.038727748  0.004439635  0.0001820986
#> [6,]  0.120801930  0.117529813  0.0681847746
#> 
#> $D
#> [1] -0.09576982 -0.07012156 -0.02410526 -0.05335663 -0.04553989 -0.09030673
#> 
```
