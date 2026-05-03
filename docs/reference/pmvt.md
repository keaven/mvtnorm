# Multivariate t Distribution

Computes the the distribution function of the multivariate t
distribution for arbitrary limits, degrees of freedom and correlation
matrices based on algorithms by Genz and Bretz.

## Usage

``` r
pmvt(lower=-Inf, upper=Inf, delta=rep(0, length(lower)),
     df=1, corr=NULL, sigma=NULL, algorithm = GenzBretz(),
     type = c("Kshirsagar", "shifted"), keepAttr=TRUE, seed = NULL, ...)
```

## Arguments

- lower:

  the vector of lower limits of length n.

- upper:

  the vector of upper limits of length n.

- delta:

  the vector of noncentrality parameters of length n, for
  `type = "shifted"` delta specifies the mode.

- df:

  degree of freedom as integer. Normal probabilities are computed for
  `df=0`.

- corr:

  the correlation matrix of dimension n.

- sigma:

  the scale matrix of dimension n. Either `corr` or `sigma` can be
  specified. If `sigma` is given, the problem is standardized
  internally. If `corr` is given, it is assumed that appropriate
  standardization was performed by the user. If neither `corr` nor
  `sigma` is given, the identity matrix is used for `sigma`.

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
  formula (1.1) in\|mvtnorm::Kotz+Nadarajah:2004\|. It is a location
  shifted version of the central t-distribution. This noncentral
  multivariate t distribution appears for example as the Bayesian
  posterior distribution for the regression coefficients in a linear
  regression. In the central case both types coincide.

- keepAttr:

  [`logical`](https://rdrr.io/r/base/logical.html) indicating if
  [`attributes`](https://rdrr.io/r/base/attributes.html) such as `error`
  and `msg` should be attached to the return value. The default, `TRUE`
  is back compatible.

- seed:

  an object specifying if and how the random number generator should be
  initialized, see [`simulate`](https://rdrr.io/r/stats/simulate.html).

- ...:

  additional parameters (currently given to `GenzBretz` for backward
  compatibility issues).

## Details

This function involves the computation of central and noncentral
multivariate t-probabilities with arbitrary correlation matrices. It
involves both the computation of singular and nonsingular probabilities.
The methodology (for default `algorithm = GenzBretz()`) is based on
randomized quasi Monte Carlo methods and described in
mvtnorm::numerical-:1999,Genz_Bretz_2002.  
Because of the randomization, the result for this algorithm (slightly)
depends on [`.Random.seed`](https://rdrr.io/r/base/Random.html).

For 2- and 3-dimensional problems one can also use the
[`TVPACK`](https://keaven.github.io/mvtnorm/reference/algorithms.md)
routines described by mvtnorm::Genz:2004, which only handles
semi-infinite integration regions (and for `type = "Kshirsagar"` only
central problems).

For `type = "Kshirsagar"` and a given correlation matrix `corr`, for
short \\A\\, say, (which has to be positive semi-definite) and degrees
of freedom \\\nu\\ the following values are numerically evaluated

\$\$I = 2^{1-\nu/2} / \Gamma(\nu/2) \int_0^\infty s^{\nu-1} \exp(-s^2/2)
\Phi(s \cdot lower/\sqrt{\nu} - \delta, s \cdot upper/\sqrt{\nu} -
\delta) \\ ds \$\$

where

\$\$\Phi(a,b) = (det(A)(2\pi)^m)^{-1/2} \int_a^b \exp(-x^\prime Ax/2) \\
dx\$\$

is the multivariate normal distribution and \\m\\ is the number of rows
of \\A\\.

For `type = "shifted"`, a positive definite symmetric matrix \\S\\
(which might be the correlation or the scale matrix), mode (vector)
\\\delta\\ and degrees of freedom \\\nu\\ the following integral is
evaluated:

\$\$c\int\_{lower_1}^{upper_1}...\int\_{lower_m}^{upper_m}
(1+(x-\delta)'S^{-1}(x-\delta)/\nu)^{-(\nu+m)/2}\\ dx_1 ... dx_m, \$\$

where

\$\$c = \Gamma((\nu+m)/2)/((\pi \nu)^{m/2}\Gamma(\nu/2)\|S\|^{1/2}),\$\$

and \\m\\ is the number of rows of \\S\\.

Note that both `-Inf` and `+Inf` may be specified in the lower and upper
integral limits in order to compute one-sided probabilities.

Univariate problems are passed to
[`pt`](https://rdrr.io/r/stats/TDist.html). If `df = 0`, normal
probabilities are returned.

## Value

The evaluated distribution function is returned, if `keepAttr` is true,
with attributes

- error:

  estimated absolute error and

- msg:

  status message (a [`character`](https://rdrr.io/r/base/character.html)
  string).

- algorithm:

  a [`character`](https://rdrr.io/r/base/character.html) string with
  `class(algorithm)`.

## References

\*

## See also

[`qmvt`](https://keaven.github.io/mvtnorm/reference/qmvt.md)

## Examples

``` r
n <- 5
lower <- -1
upper <- 3
df <- 4
corr <- diag(5)
corr[lower.tri(corr)] <- 0.5
delta <- rep(0, 5)
prob <- pmvt(lower=lower, upper=upper, delta=delta, df=df, corr=corr)
print(prob)
#> [1] 0.5064858
#> attr(,"error")
#> [1] 0.0003305051
#> attr(,"msg")
#> [1] "Normal Completion"

pmvt(lower=-Inf, upper=3, df = 3, sigma = 1) == pt(3, 3)
#> upper 
#>  TRUE 

# Example from R News paper (original by Edwards and Berry, 1987)

n <- c(26, 24, 20, 33, 32)
V <- diag(1/n)
df <- 130
C <- c(1,1,1,0,0,-1,0,0,1,0,0,-1,0,0,1,0,0,0,-1,-1,0,0,-1,0,0)
C <- matrix(C, ncol=5)
### scale matrix
cv <- C %*% tcrossprod(V, C)
### correlation matrix
cr <- cov2cor(cv)
delta <- rep(0,5)

myfct <- function(q, alpha) {
  lower <- rep(-q, ncol(cv))
  upper <- rep(q, ncol(cv))
  pmvt(lower=lower, upper=upper, delta=delta, df=df,
       corr=cr, abseps=0.0001) - alpha
}

### uniroot for this simple problem
round(uniroot(myfct, lower=1, upper=5, alpha=0.95)$root, 3)
#> [1] 2.561

# compare pmvt and pmvnorm for large df:

a <- pmvnorm(lower=-Inf, upper=1, mean=rep(0, 5), corr=diag(5))
b <- pmvt(lower=-Inf, upper=1, delta=rep(0, 5), df=300,
          corr=diag(5))
a
#> [1] 0.4215702
#> attr(,"error")
#> [1] 0
#> attr(,"msg")
#> [1] "Normal Completion"
b
#> [1] 0.421141
#> attr(,"error")
#> [1] 1.635956e-06
#> attr(,"msg")
#> [1] "Normal Completion"

stopifnot(round(a, 2) == round(b, 2))

# correlation and scale matrix

a <- pmvt(lower=-Inf, upper=2, delta=rep(0,5), df=3,
          sigma = diag(5)*2)
b <- pmvt(lower=-Inf, upper=2/sqrt(2), delta=rep(0,5),
          df=3, corr=diag(5))
attributes(a) <- NULL
attributes(b) <- NULL
a
#> [1] 0.5654061
b
#> [1] 0.5654025
stopifnot(all.equal(round(a,3) , round(b, 3)))

a <- pmvt(0, 1,df=10)
attributes(a) <- NULL
b <- pt(1, df=10) - pt(0, df=10)
stopifnot(all.equal(round(a,10) , round(b, 10)))
```
