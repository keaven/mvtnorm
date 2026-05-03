# Multivariate Normal Distribution

Computes the distribution function of the multivariate normal
distribution for arbitrary limits and correlation matrices.

## Usage

``` r
pmvnorm(lower=-Inf, upper=Inf, mean=rep(0, length(lower)),
           corr=NULL, sigma=NULL, algorithm = GenzBretz(), keepAttr=TRUE, 
           seed = NULL, ...)
```

## Arguments

- lower:

  the vector of lower limits of length n.

- upper:

  the vector of upper limits of length n.

- mean:

  the mean vector of length n.

- corr:

  the correlation matrix of dimension n.

- sigma:

  the covariance matrix of dimension n less than 1000. Either `corr` or
  `sigma` can be specified. If `sigma` is given, the problem is
  standardized internally. If `corr` is given, it is assumed that
  appropriate standardization was performed by the user. If neither
  `corr` nor `sigma` is given, the identity matrix is used for `sigma`.

- algorithm:

  an object of class
  [`GenzBretz`](https://keaven.github.io/mvtnorm/reference/algorithms.md),
  [`Miwa`](https://keaven.github.io/mvtnorm/reference/algorithms.md) or
  [`TVPACK`](https://keaven.github.io/mvtnorm/reference/algorithms.md)
  specifying both the algorithm to be used as well as the associated
  hyper parameters.

- keepAttr:

  [`logical`](https://rdrr.io/r/base/logical.html) indicating if
  [`attributes`](https://rdrr.io/r/base/attributes.html) such as `error`
  and `msg` should be attached to the return value. The default, `TRUE`
  is back compatible.

- seed:

  an object specifying if and how the random number generator should be
  initialized, see [`simulate`](https://rdrr.io/r/stats/simulate.html).

- ...:

  additional parameters (currently given to
  [`GenzBretz`](https://keaven.github.io/mvtnorm/reference/algorithms.md)
  for backward compatibility issues).

## Details

This program involves the computation of multivariate normal
probabilities with arbitrary correlation matrices. It involves both the
computation of singular and nonsingular probabilities. The implemented
methodology is described in mvtnorm::numerical-:1992 and
mvtnorm::comparison:1993 for algorithm `GenzBretz`, in
mvtnorm::Miwa+Hayter+Kuriki:2003 for algorithm `Miwa`, useful up to
dimension 20, and mvtnorm::Genz:2004 for the `TVPACK` algorithm, which
covers 2- and 3-dimensional problems for semi-infinite integration
regions.

Note the default algorithm GenzBretz is randomized and hence slightly
depends on [`.Random.seed`](https://rdrr.io/r/base/Random.html) and that
both `-Inf` and `+Inf` may be specified in `lower` and `upper`. For more
details see
[`pmvt`](https://keaven.github.io/mvtnorm/reference/pmvt.md).

The multivariate normal case is treated as a special case of
[`pmvt`](https://keaven.github.io/mvtnorm/reference/pmvt.md) with `df=0`
and univariate problems are passed to
[`pnorm`](https://rdrr.io/r/stats/Normal.html).

The multivariate normal density and random deviates are available using
[`dmvnorm`](https://keaven.github.io/mvtnorm/reference/Mvnorm.md) and
[`rmvnorm`](https://keaven.github.io/mvtnorm/reference/Mvnorm.md).

`pmvnorm` is based on original implementations by Alan Genz, Frank
Bretz, and Tetsuhisa Miwa developed for computing accurate
approximations to the normal integral. Users interested in computing
log-likelihoods involving such normal probabilities should consider
function
[`lpmvnorm`](https://keaven.github.io/mvtnorm/reference/lpmvnorm.md),
which is more flexible and efficient for this task and comes with the
ability to evaluate score functions.

An overview is available from mvtnorm::Genz_Bretz_2009.

## Value

The evaluated distribution function is returned, if `keepAttr` is true,
with attributes

- error:

  estimated absolute error

- msg:

  status message(s).

- algorithm:

  a [`character`](https://rdrr.io/r/base/character.html) string with
  `class(algorithm)`.

## References

\*

## See also

[`qmvnorm`](https://keaven.github.io/mvtnorm/reference/qmvnorm.md) for
quantiles and
[`lpmvnorm`](https://keaven.github.io/mvtnorm/reference/lpmvnorm.md) for
log-likelihoods.

## Examples

``` r
n <- 5
mean <- rep(0, 5)
lower <- rep(-1, 5)
upper <- rep(3, 5)
corr <- diag(5)
corr[lower.tri(corr)] <- 0.5
corr[upper.tri(corr)] <- 0.5
prob <- pmvnorm(lower, upper, mean, corr)
print(prob)
#> [1] 0.5801087
#> attr(,"error")
#> [1] 0.0003069052
#> attr(,"msg")
#> [1] "Normal Completion"

stopifnot(pmvnorm(lower=-Inf, upper=3, mean=0, sigma=1) == pnorm(3))

a <- pmvnorm(lower=-Inf,upper=c(.3,.5),mean=c(2,4),diag(2))

stopifnot(round(a,16) == round(prod(pnorm(c(.3,.5),c(2,4))),16))

a <- pmvnorm(lower=-Inf,upper=c(.3,.5,1),mean=c(2,4,1),diag(3))

stopifnot(round(a,16) == round(prod(pnorm(c(.3,.5,1),c(2,4,1))),16))

# Example from R News paper (original by Genz, 1992):

m <- 3
sigma <- diag(3)
sigma[2,1] <- 3/5
sigma[3,1] <- 1/3
sigma[3,2] <- 11/15
pmvnorm(lower=rep(-Inf, m), upper=c(1,4,2), mean=rep(0, m), corr=sigma)
#> [1] 0.8279848
#> attr(,"error")
#> [1] 5.20711e-07
#> attr(,"msg")
#> [1] "Normal Completion"

# Correlation and Covariance

a <- pmvnorm(lower=-Inf, upper=c(2,2), sigma = diag(2)*2)
b <- pmvnorm(lower=-Inf, upper=c(2,2)/sqrt(2), corr=diag(2))
stopifnot(all.equal(round(a,5) , round(b, 5)))
```
