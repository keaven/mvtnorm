# mvtnorm 1.3-7.9000

Development version for the `Deterministic-QMC-algorithm` branch. The
package's historical CRAN release notes remain in `inst/NEWS.Rd`.

### New

- Added `DeterministicQMC()` as a reproducible alternative to randomized
  `GenzBretz()` integration.
- Added deterministic C backends for multivariate normal probabilities:
  `backend = "c"`, `"sobol"`, `"hybrid"`, and `"genzsobol"`.
- Added `backend = "auto"` to select a deterministic policy without user
  intervention. Supported normal problems use cross-checked Genz-Sobol;
  unsupported cases use the deterministic FORTRAN path.
- Added diagnostic tolerance handling for deterministic backends. Positive
  `abseps` or `releps` requests are checked against deterministic stability
  diagnostics and return `inform = 1` if the requested diagnostic target is not
  met.
- Added deterministic QMC regression tests.
- Added pkgdown configuration and deterministic QMC articles:
  - `vignettes/articles/dqmc-methods.Rmd`
  - `vignettes/articles/dqmc-benchmark.Rmd`
  - matching Sweave vignettes under `vignettes/`
- Added generated deterministic QMC vignette artifacts under `inst/doc/`.

### Notes

- The deterministic `error` attribute is a reproducible stability diagnostic,
  not a Monte Carlo confidence interval and not a rigorous deterministic error
  bound.
- The current Sobol implementation uses an embedded plain C direction-number
  table rather than modern high-dimensional Joe-Kuo direction files.
- The `metal` backend currently probes for a macOS Metal device but falls back
  to C/FORTRAN because no Metal evaluator is implemented yet.
- `backend = "auto"` currently uses fixed deterministic sequence offsets. A
  future implementation should consider stronger deterministic digital shifts
  and more mature adaptive refinement.
