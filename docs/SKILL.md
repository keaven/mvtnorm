# mvtnorm Deterministic QMC Branch

Use this project guide when editing this repository. The current branch
adds an experimental deterministic QMC backend to `mvtnorm` while
keeping the original FORTRAN/C algorithms as reference implementations.

## Current Branch

- Branch: `Deterministic-QMC-algorithm`
- Upstream remote: `https://github.com/cran/mvtnorm.git`
- Primary new user API:
  [`DeterministicQMC()`](https://keaven.github.io/mvtnorm/reference/algorithms.md)
- Main R implementation file: `R/mvt.R`
- Main C implementation file: `src/dqmc.c`
- Native registration files:
  - `src/mvtnorm-init.c`
  - `src/mvtnorm.h`
  - `src/C_FORTRAN_interface.c`

## Deterministic QMC Backends

- `backend = "c"`: deterministic Halton C backend for normal
  probabilities.
- `backend = "sobol"`: deterministic Sobol C backend.
- `backend = "hybrid"`: Sobol plus deterministic tail-spaced warp.
- `backend = "genzsobol"`: Sobol plus variable reordering and embedded
  stability checks.
- `backend = "auto"`: supported normal probabilities use repeated
  Genz-Sobol checks at fixed sequence offsets; unsupported cases use
  deterministic FORTRAN.
- `backend = "metal"`: currently probes for Metal availability but falls
  back because no Metal evaluator is implemented.

## Documentation To Keep In Sync

- Help file: `man/algorithms.Rd`
- Methodology article: `vignettes/articles/dqmc-methods.Rmd`
- Benchmark article: `vignettes/articles/dqmc-benchmark.Rmd`
- Sweave vignettes:
  - `vignettes/DeterministicQMC-methods.Rnw`
  - `vignettes/DeterministicQMC-benchmark.Rnw`
- Generated installed docs under `inst/doc/`
- pkgdown site under `docs/`
- `MD5`

## Accuracy Language

Use careful wording for deterministic error control. `abseps` and
`releps` request reproducible stability diagnostics for
[`DeterministicQMC()`](https://keaven.github.io/mvtnorm/reference/algorithms.md).
They are not stochastic confidence intervals and not rigorous
deterministic error bounds.

## Useful Commands

Install locally into the temporary test library:

``` sh
R CMD INSTALL -l /tmp/mvtnorm-lib .
```

Run the deterministic QMC regression test:

``` sh
R_LIBS=/tmp/mvtnorm-lib Rscript tests/regtest-DeterministicQMC.R
```

Build deterministic QMC Sweave vignettes:

``` sh
R_LIBS=/tmp/mvtnorm-lib Rscript -e 'tools::buildVignette("vignettes/DeterministicQMC-benchmark.Rnw", dir = ".", tangle = TRUE); tools::buildVignette("vignettes/DeterministicQMC-methods.Rnw", dir = ".", tangle = TRUE)'
```

Build pkgdown:

``` sh
Rscript -e 'pkgdown::build_site(new_process = FALSE)'
```

Run a package check without rebuilding vignettes:

``` sh
mkdir -p /tmp/mvtnorm-check
R CMD check --no-manual --no-build-vignettes -o /tmp/mvtnorm-check .
```

Update and verify `MD5` after source/doc changes:

``` sh
Rscript -e 'paths <- substring(readLines("MD5"), 35L); sums <- tools::md5sum(paths); writeLines(sprintf("%s *%s", unname(sums), paths), "MD5")'
Rscript -e 'x <- readLines("MD5"); paths <- substring(x, 35L); sums <- tools::md5sum(paths); stopifnot(identical(unname(sums), substr(x, 1L, 32L))); cat("MD5 OK\n")'
```

Clean native build artifacts before committing:

``` sh
rm -f src/C_FORTRAN_interface.o src/dqmc.o src/lpmvnorm.o src/ltMatrices.o src/miwa.o src/mvt.o src/mvtnorm-init.o src/mvtnorm.so src/tvpack.o src/symbols.rds
```

## Known Check Noise

The current branch has known `R CMD check` warnings from Rd bibliography
macros such as `\bibcitet` and `\bibcitep`, plus environment-specific
compiler header warnings. Hidden local files such as `.git` and
`.Rproj.user` can produce a NOTE when checking the source directory
directly.
