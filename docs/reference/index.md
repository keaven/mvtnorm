# Package index

## Multivariate normal functions

Core package overview, distribution functions, likelihoods, and
normal-model interfaces.

- [`mvtnorm-package`](https://keaven.github.io/mvtnorm/reference/mvtnorm-package.md)
  [`mvtnorm`](https://keaven.github.io/mvtnorm/reference/mvtnorm-package.md)
  : Multivariate Normal and t Distributions
- [`pmvnorm()`](https://keaven.github.io/mvtnorm/reference/pmvnorm.md) :
  Multivariate Normal Distribution
- [`dmvnorm()`](https://keaven.github.io/mvtnorm/reference/Mvnorm.md)
  [`rmvnorm()`](https://keaven.github.io/mvtnorm/reference/Mvnorm.md) :
  Multivariate Normal Density and Random Deviates
- [`qmvnorm()`](https://keaven.github.io/mvtnorm/reference/qmvnorm.md) :
  Quantiles of the Multivariate Normal Distribution
- [`lpmvnorm()`](https://keaven.github.io/mvtnorm/reference/lpmvnorm.md)
  [`slpmvnorm()`](https://keaven.github.io/mvtnorm/reference/lpmvnorm.md)
  [`ldmvnorm()`](https://keaven.github.io/mvtnorm/reference/lpmvnorm.md)
  [`sldmvnorm()`](https://keaven.github.io/mvtnorm/reference/lpmvnorm.md)
  [`ldpmvnorm()`](https://keaven.github.io/mvtnorm/reference/lpmvnorm.md)
  [`sldpmvnorm()`](https://keaven.github.io/mvtnorm/reference/lpmvnorm.md)
  : Multivariate Normal Log-likelihood and Score Functions
- [`lpRR()`](https://keaven.github.io/mvtnorm/reference/lpRR.md)
  [`slpRR()`](https://keaven.github.io/mvtnorm/reference/lpRR.md) :
  Multivariate Normal Log-likelihood and Score Functions for Reduced
  Rank Covariances
- [`mvnorm()`](https://keaven.github.io/mvtnorm/reference/interface.md)
  [`aperm(`*`<mvnorm>`*`)`](https://keaven.github.io/mvtnorm/reference/interface.md)
  [`margDist()`](https://keaven.github.io/mvtnorm/reference/interface.md)
  [`condDist()`](https://keaven.github.io/mvtnorm/reference/interface.md)
  [`simulate(`*`<mvnorm>`*`)`](https://keaven.github.io/mvtnorm/reference/interface.md)
  [`logLik(`*`<mvnorm>`*`)`](https://keaven.github.io/mvtnorm/reference/interface.md)
  [`lLgrad(`*`<mvnorm>`*`)`](https://keaven.github.io/mvtnorm/reference/interface.md)
  : User Interface to Multiple Multivariate Normal Distributions
- [`marg_mvnorm()`](https://keaven.github.io/mvtnorm/reference/margcond.md)
  [`cond_mvnorm()`](https://keaven.github.io/mvtnorm/reference/margcond.md)
  : Marginal and Conditional Multivariate Normal Distributions

## t-distribution functions

Probability, density, random generation, and quantile functions for
multivariate t distributions.

- [`pmvt()`](https://keaven.github.io/mvtnorm/reference/pmvt.md) :
  Multivariate t Distribution
- [`rmvt()`](https://keaven.github.io/mvtnorm/reference/Mvt.md)
  [`dmvt()`](https://keaven.github.io/mvtnorm/reference/Mvt.md) : The
  Multivariate t Distribution
- [`qmvt()`](https://keaven.github.io/mvtnorm/reference/qmvt.md) :
  Quantiles of the Multivariate t Distribution

## Algorithm and hyper parameter selection

Algorithm constructors and integration-control parameters.

- [`GenzBretz()`](https://keaven.github.io/mvtnorm/reference/algorithms.md)
  [`DeterministicQMC()`](https://keaven.github.io/mvtnorm/reference/algorithms.md)
  [`Miwa()`](https://keaven.github.io/mvtnorm/reference/algorithms.md)
  [`TVPACK()`](https://keaven.github.io/mvtnorm/reference/algorithms.md)
  : Choice of Algorithm and Hyper Parameters

## Matrix functions

Lower-triangular and symmetric matrix classes, conversions, and linear
algebra helpers.

- [`ltMatrices()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`syMatrices()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`as.array(`*`<ltMatrices>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`as.array(`*`<syMatrices>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`diagonals(`*`<ltMatrices>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`diagonals(`*`<syMatrices>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`diagonals(`*`<matrix>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`diagonals(`*`<integer>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`` `diagonals<-`() ``](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`solve(`*`<ltMatrices>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`chol(`*`<syMatrices>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`aperm(`*`<chol>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`aperm(`*`<invchol>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`aperm(`*`<ltMatrices>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`aperm(`*`<syMatrices>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`deperma()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`Mult(`*`<ltMatrices>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`Mult(`*`<syMatrices>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`Tcrossprod()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`Crossprod()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`tcrossprod(`*`<ltMatrices>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`tcrossprod(`*`<syMatrices>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`crossprod(`*`<ltMatrices>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`crossprod(`*`<syMatrices>`*`)`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`logdet()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`Lower_tri()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`is.ltMatrices()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`is.syMatrices()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`as.ltMatrices()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`as.syMatrices()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`is.chol()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`is.invchol()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`as.chol()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`as.invchol()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`chol2cov()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`invchol2chol()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`chol2invchol()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`invchol2cov()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`invchol2pre()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`chol2pre()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`Dchol()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`invcholD()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`chol2cor()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`invchol2cor()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`chol2pc()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`invchol2pc()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`vectrick()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`standardize()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  [`destandardize()`](https://keaven.github.io/mvtnorm/reference/ltMatrices.md)
  : Multiple Lower Triangular or Symmetric Matrices
