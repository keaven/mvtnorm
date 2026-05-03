# mvtnorm 1.3-7-9001

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
- Added development tooling for the fork:
  - GitHub Actions `R CMD check` across Linux, macOS, Windows, and multiple
    R versions.
  - A pkgdown build workflow with GitHub Pages deployment on default branches.
  - A coverage workflow using `covr`.
  - An incremental `testthat` edition 3 test suite for new deterministic QMC
    behavior while retaining the existing regression tests.
- Added Sobol, quasi-Monte Carlo, and Genz integration references for the new
  deterministic QMC backends.

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
