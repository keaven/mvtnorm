# Choice of Algorithm and Hyper Parameters

Choose between three algorithms for evaluating normal (and t-)
distributions and define hyper parameters.

## Usage

``` r
GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0)
DeterministicQMC(maxpts = 25000, abseps = 0.001, releps = 0,
                 backend = c("auto", "c", "sobol", "hybrid",
                             "genzsobol", "fortran", "cpu", "metal"),
                 start = 0L)
Miwa(steps = 128, checkCorr = TRUE, maxval = 1e3)
TVPACK(abseps = 1e-6)
```

## Arguments

- maxpts:

  maximum number of function values as integer. The internal FORTRAN
  code always uses a minimum number depending on the dimension. (for
  example 752 for three-dimensional problems).

- abseps:

  absolute error tolerance. For `DeterministicQMC`, this is a
  deterministic stability-diagnostic target, not a stochastic confidence
  bound. For `TVPACK`, it is only used for dimension 3.

- releps:

  relative error tolerance as double. For `DeterministicQMC`, this is
  interpreted as a deterministic stability-diagnostic target.

- backend:

  computational backend for `DeterministicQMC`. `"auto"` selects a
  deterministic policy without user intervention: supported normal
  probabilities use a cross-checked Genz-Sobol path, while
  higher-dimensional normal probabilities and t-probabilities use the
  deterministic FORTRAN implementation. `"c"` uses a deterministic
  Halton rule for normal probabilities, and `"sobol"` uses a
  deterministic Sobol digital net for normal probabilities. `"hybrid"`
  is an experimental Sobol backend with a deterministic tail-spaced
  change of variables. `"genzsobol"` is an experimental Sobol backend
  with marginal-width variable reordering and embedded powers-of-two
  stability checks. `"cpu"` is an alias for `"c"`. The `"metal"` option
  probes for a default macOS Metal device and currently falls back to
  C/FORTRAN.

- start:

  non-negative integer offset for the deterministic shift sequence used
  by `DeterministicQMC`.

- steps:

  number of grid points to be evaluated; cannot be larger than 4097.

- checkCorr:

  logical indicating if a check for singularity of the correlation
  matrix should be performed (once per function call to
  [`pmvt()`](https://keaven.github.io/mvtnorm/reference/pmvt.md) or
  [`pmvnorm()`](https://keaven.github.io/mvtnorm/reference/pmvnorm.md)).

- maxval:

  replacement for `Inf` when non-orthant probabilities involving `Inf`
  shall be computed.

## Details

There are three algorithms available for evaluating normal (and two
algorithms for t-) probabilities: The default is the randomized
Quasi-Monte-Carlo procedure by
mvtnorm::numerical-:1992,mvtnorm::comparison:1993 and
mvtnorm::Genz_Bretz_2002 applicable to arbitrary covariance structures
and dimensions up to 1000.

`DeterministicQMC` uses the same transformed integration problem as
`GenzBretz`, but replaces randomized integration by deterministic
quasi-Monte Carlo rules. The `"c"` backend currently implements normal
probabilities using a deterministic Halton rule. The `"sobol"` backend
implements the same plain C transformed normal integration with Sobol
direction numbers and Gray-code point generation; it currently supports
normal probabilities up to 33 dimensions and falls back to the Halton C
backend for higher dimensions. The experimental `"hybrid"` backend uses
the Sobol point set with a deterministic tail-spaced change of variables
and the corresponding Jacobian adjustment, to test whether sparse tail
spacing improves fixed-budget accuracy. The experimental `"genzsobol"`
backend adds marginal-width variable reordering and embedded
powers-of-two stability checks using `abseps` and `releps`. With
`backend = "auto"`, supported normal probabilities are evaluated by
repeated deterministic Genz-Sobol checks at fixed sequence offsets. When
`abseps` or `releps` is requested, the individual checks use a tighter
internal target before the combined diagnostic is compared with the
user's requested tolerance. The returned value is the mean of these
checks and the reported error is the larger of the embedded stability
errors and the between-check spread. If this deterministic diagnostic
exceeds `abseps` or `releps`, the result is returned with `inform = 1`;
no stochastic fallback is used. The same diagnostic target is applied to
the reported deterministic block-stability error for the other C
backends. **Warning:** `abseps` and `releps` for `DeterministicQMC` are
reproducible stability checks. They should not be interpreted as Monte
Carlo confidence statements or rigorous deterministic error bounds. The
`"fortran"` backend uses the original FORTRAN implementation with
deterministic shifts. This makes repeated evaluations reproducible
without setting R's random-number seed. The reported error is a
deterministic convergence diagnostic rather than a randomized Monte
Carlo confidence statement.

For normal probabilities, smaller dimensions (up to 20) and non-singular
covariance matrices, the algorithm by mvtnorm::Miwa+Hayter+Kuriki:2003
can be used as well. This algorithm can compute orthant probabilities
(`lower` being `-Inf` or `upper` equal to `Inf`). Non-orthant
probabilities are computed from the corresponding orthant probabilities,
however, infinite limits are replaced by `maxval` along with a warning.

For two- and three-dimensional problems and semi-infinite integration
region, `TVPACK` implements an interface to the methods described by
mvtnorm::Genz:2004.

## Value

An object of class `"GenzBretz"`, `"DeterministicQMC"`, `"Miwa"`, or
`"TVPACK"` defining hyper parameters.

## References

\*
