# Multivariate Normal and t Distributions

Computes multivariate normal and t probabilities, quantiles, random
deviates, and densities. Log-likelihoods for multivariate Gaussian
models and Gaussian copulae parameterised by Cholesky factors of
covariance or precision matrices are implemented for interval-censored
and exact data, or a mix thereof. Score functions for these
log-likelihoods are available. A class representing multiple lower
triangular matrices and corresponding methods are part of this package.

## Details

Package mvtnorm provides functionality for dealing with multivariate
normal and t-distributions. The package interfaces FORTRAN and `C` code
for evaluating multivariate normal probabilities written by Alan Genz
and Tetsuhisa Miwa. Functions
[`pmvnorm`](http://mvtnorm.R-forge.R-project.org/reference/pmvnorm.md),
[`pmvt`](http://mvtnorm.R-forge.R-project.org/reference/pmvt.md),
[`qmvnorm`](http://mvtnorm.R-forge.R-project.org/reference/qmvnorm.md),
and [`qmvt`](http://mvtnorm.R-forge.R-project.org/reference/qmvt.md)
return normal and t probabilities or corresponding quantiles computed by
these original implementations. Users interested in the computation of
such probabilities or quantiles, for example for multiple testing
purposes, should use this functionality.

When the multivariate normal log-likelihood function, defined by the
log-probability in the discrete or interval-censored case or by the
log-density for exact real observations, or a mix thereof, shall be
computed, functions
[`lpmvnorm`](http://mvtnorm.R-forge.R-project.org/reference/lpmvnorm.md),
[`ldmvnorm`](http://mvtnorm.R-forge.R-project.org/reference/lpmvnorm.md),
and
[`ldpmvnorm`](http://mvtnorm.R-forge.R-project.org/reference/lpmvnorm.md)
are better suited. They rely on an independent implementation of Genz'
algorithm (for log-probabilities), can be customised (different
quasi-Monte Carlo schemes), and are a bit faster. Most importantly, the
corresponding score functions are available through functions
[`slpmvnorm`](http://mvtnorm.R-forge.R-project.org/reference/lpmvnorm.md),
[`sldmvnorm`](http://mvtnorm.R-forge.R-project.org/reference/lpmvnorm.md),
or
[`sldpmvnorm`](http://mvtnorm.R-forge.R-project.org/reference/lpmvnorm.md),
which help to speed-up parameter estimation considerably. Users
interested in this functionality should consult the `lmvnorm_src`
package vignette.

## See also

[`vignette("lmvnorm_src", package = "mvtnorm")`](https://cran.rstudio.com/web/packages/mvtnorm/vignettes/lmvnorm_src.pdf)
