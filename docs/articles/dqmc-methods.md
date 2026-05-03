# Deterministic QMC Methodology

## The computational problem

[`pmvnorm()`](https://keaven.github.io/mvtnorm/reference/pmvnorm.md) and
[`pmvt()`](https://keaven.github.io/mvtnorm/reference/pmvt.md) evaluate
probabilities over rectangular regions for multivariate normal and
Student t distributions. After standardisation, the problem is an
integral over a unit cube or hyperrectangle whose dimension is the
dimension of the distribution. Tensor-product quadrature becomes
unattractive quickly, so `mvtnorm` uses specialised transformations and
simulation or quasi-Monte Carlo rules.

The default
[`GenzBretz()`](https://keaven.github.io/mvtnorm/reference/algorithms.md)
algorithm is randomized quasi-Monte Carlo. Its randomization is a
feature: it makes repeated randomized shifts possible and therefore
gives a practical stochastic error assessment while retaining the speed
advantages of low-discrepancy integration points.

The new
[`DeterministicQMC()`](https://keaven.github.io/mvtnorm/reference/algorithms.md)
option keeps the same broad Genz-style transformed integration
viewpoint, but removes dependence on R’s random-number generator. The
plain C backends implement normal probabilities with a Cholesky
factorization and sequential conditional integration. `backend = "c"`
uses a deterministic Halton rule; `backend = "sobol"` uses a
deterministic Sobol digital net; `backend = "hybrid"` tests a Sobol rule
with a deterministic tail-spaced change of variables;
`backend = "genzsobol"` adds marginal-width variable reordering and
embedded stability checks. The default deterministic choice,
`backend = "auto"`, now selects a cross-checked Genz-Sobol policy for
supported normal problems and keeps the deterministic FORTRAN path for
higher-dimensional or Student t problems. The original FORTRAN code
remains as a reference implementation and fallback, including for
Student t probabilities.

## Method comparison

| Algorithm | Main advantages | Main limitations |
|----|----|----|
| [`GenzBretz()`](https://keaven.github.io/mvtnorm/reference/algorithms.md) | Mature default; fast randomized QMC; supports normal and t probabilities; randomized error assessment. | Repeated calls can differ unless the seed is controlled; stochastic noise can be awkward in deterministic optimization and regression tests. |
| [`Miwa()`](https://keaven.github.io/mvtnorm/reference/algorithms.md) | Deterministic; often very accurate for low-dimensional normal probabilities; independent of RNG state. | Normal probabilities only; dimension limited to 20; grid cost grows quickly. |
| Fixed-budget original QMC | Uses the established randomized Genz-Bretz path with early stopping disabled, making the point budget explicit. | Still randomized; error assessment is tied to randomized replication; runtime depends directly on the chosen budget. |
| [`DeterministicQMC()`](https://keaven.github.io/mvtnorm/reference/algorithms.md) | Deterministic output without managing seeds; plain C Halton, Sobol, hybrid, and Genz-Sobol normal backends; FORTRAN reference/fallback remains available; uses `maxpts` budget control. | The C backends currently handle normal probabilities only; their error values are deterministic convergence diagnostics rather than randomized confidence statements. |

## What is new and useful

The new backend is not intended to discard the FORTRAN implementation.
It adds a second deterministic route through the existing `mvtnorm` API:

- `DeterministicQMC(backend = "c")` calls a plain C evaluator for
  multivariate normal probabilities using a deterministic Halton rule.
- `DeterministicQMC(backend = "sobol")` calls the same plain C
  transformed normal evaluator with Sobol direction numbers and
  Gray-code point generation.
- `DeterministicQMC(backend = "hybrid")` uses the Sobol generator with a
  deterministic tail-spaced change of variables and Jacobian correction.
  This is experimental and intended to test whether sparse tail spacing
  improves fixed-budget accuracy.
- `DeterministicQMC(backend = "genzsobol")` adds two algorithmic pieces
  from the proposed fuller backend: variables are reordered by marginal
  probability width before the Cholesky factorization, and estimates are
  checked at embedded powers of two against `abseps`/`releps`.
- `DeterministicQMC(backend = "fortran")` keeps the historical FORTRAN
  code as a reference/fallback path with deterministic shifts.
- `DeterministicQMC(backend = "auto")` is the best-practice
  deterministic selector: normal probabilities supported by the C Sobol
  code use Genz-Sobol with fixed deterministic cross-checks, while
  unsupported normal problems and t probabilities use deterministic
  FORTRAN rather than a stochastic fallback.
- `DeterministicQMC(backend = "metal")` probes for a macOS Metal device,
  but currently falls back because a Metal evaluator has not yet been
  implemented.

The main advantage over the previous deterministic option,
[`Miwa()`](https://keaven.github.io/mvtnorm/reference/algorithms.md), is
not a claim of universal accuracy. The advantage is architectural:
deterministic QMC scales by a point budget (`maxpts`) rather than by a
low-dimensional grid. That makes it more natural for deterministic
optimization, reproducible simulation studies, and regression tests
where exact repeatability is important.

For paper citations of the new deterministic Sobol and Genz-Sobol
backends, a good citation trail is:

- Sobol (1967) for the Sobol digital sequence;
- Caflisch (1998) or Glasserman (2004, Chapter 6) for quasi-Monte Carlo
  background, convergence motivation, and practical limitations;
- Genz (1992), Genz and Bretz (2002, 2009), and the existing `mvtnorm`
  references for the transformed multivariate normal and Student t
  integration algorithms.

The corresponding BibTeX entries are kept in `inst/REFERENCES.bib`.

``` r
GenzBretz()
```

    ## $maxpts
    ## [1] 25000
    ## 
    ## $abseps
    ## [1] 0.001
    ## 
    ## $releps
    ## [1] 0
    ## 
    ## attr(,"class")
    ## [1] "GenzBretz"

``` r
DeterministicQMC(backend = "auto")
```

    ## $maxpts
    ## [1] 25000
    ## 
    ## $abseps
    ## [1] 0.001
    ## 
    ## $releps
    ## [1] 0
    ## 
    ## $backend
    ## [1] "auto"
    ## 
    ## $start
    ## [1] 0
    ## 
    ## attr(,"class")
    ## [1] "DeterministicQMC" "GenzBretz"

``` r
DeterministicQMC(backend = "c")
```

    ## $maxpts
    ## [1] 25000
    ## 
    ## $abseps
    ## [1] 0.001
    ## 
    ## $releps
    ## [1] 0
    ## 
    ## $backend
    ## [1] "c"
    ## 
    ## $start
    ## [1] 0
    ## 
    ## attr(,"class")
    ## [1] "DeterministicQMC" "GenzBretz"

``` r
DeterministicQMC(backend = "sobol")
```

    ## $maxpts
    ## [1] 25000
    ## 
    ## $abseps
    ## [1] 0.001
    ## 
    ## $releps
    ## [1] 0
    ## 
    ## $backend
    ## [1] "sobol"
    ## 
    ## $start
    ## [1] 0
    ## 
    ## attr(,"class")
    ## [1] "DeterministicQMC" "GenzBretz"

``` r
DeterministicQMC(backend = "hybrid")
```

    ## $maxpts
    ## [1] 25000
    ## 
    ## $abseps
    ## [1] 0.001
    ## 
    ## $releps
    ## [1] 0
    ## 
    ## $backend
    ## [1] "hybrid"
    ## 
    ## $start
    ## [1] 0
    ## 
    ## attr(,"class")
    ## [1] "DeterministicQMC" "GenzBretz"

``` r
DeterministicQMC(backend = "genzsobol")
```

    ## $maxpts
    ## [1] 25000
    ## 
    ## $abseps
    ## [1] 0.001
    ## 
    ## $releps
    ## [1] 0
    ## 
    ## $backend
    ## [1] "genzsobol"
    ## 
    ## $start
    ## [1] 0
    ## 
    ## attr(,"class")
    ## [1] "DeterministicQMC" "GenzBretz"

``` r
DeterministicQMC(backend = "fortran")
```

    ## $maxpts
    ## [1] 25000
    ## 
    ## $abseps
    ## [1] 0.001
    ## 
    ## $releps
    ## [1] 0
    ## 
    ## $backend
    ## [1] "fortran"
    ## 
    ## $start
    ## [1] 0
    ## 
    ## attr(,"class")
    ## [1] "DeterministicQMC" "GenzBretz"

``` r
Miwa()
```

    ## $steps
    ## [1] 128
    ## 
    ## $checkCorr
    ## [1] TRUE
    ## 
    ## $maxval
    ## [1] 1000
    ## 
    ## attr(,"class")
    ## [1] "Miwa"

## The Genz-Sobol backends and the faster target

The first C backend deliberately used a simple deterministic Halton
rule. The new `backend = "sobol"` implements the next point-generator
step: it keeps the sequential Genz transformation and replaces Halton
points with an embedded Sobol digital net generated from direction
numbers in plain C. The experimental `backend = "hybrid"` keeps that
Sobol backbone and applies a gentle deterministic tail-spaced warp with
a Jacobian adjustment. This lets the benchmark test the nonuniform-grid
idea without replacing QMC by a tensor grid. The `backend = "genzsobol"`
option keeps the Sobol backbone, reorders variables before
factorization, and can stop early when embedded powers-of-two estimates
stabilize.

The backend is still intentionally conservative. It does not yet include
all of the refinements that make the mature randomized Genz-Bretz path
so effective, but `backend = "auto"` now implements the first automatic
policy: run Genz-Sobol at several fixed deterministic sequence offsets,
combine the estimates, and use the larger of embedded convergence error
and between-check spread as the accuracy diagnostic. When a tolerance is
requested, the repeated checks are run against a tighter internal target
before the combined diagnostic is compared with the user-requested
tolerance. The remaining speed-and-accuracy target is to add:

- reorder variables before integration so the transformed integrand has
  lower effective dimension;
- compare the Sobol backend with the experimental tail-spaced hybrid;
- replace the present fixed sequence offsets by stronger deterministic
  digital shifts or rank-1 lattice shifts while preserving exact
  reproducibility;
- adaptively double the point budget until embedded estimates stabilize;
- fall back to the deterministic FORTRAN reference path, or to a larger
  budget, when the deterministic stability diagnostic is poor.

The Sobol backend does not have the same randomized error interpretation
as
[`GenzBretz()`](https://keaven.github.io/mvtnorm/reference/algorithms.md),
because that interpretation comes from randomization. Its goal is
instead comparable numerical accuracy with deterministic output. This is
plausible because it uses the same Genz-style conditioning strategy and
a stronger low-discrepancy design than plain Halton; further
effective-dimension reduction remains a natural enhancement. The hybrid
is useful to test because it preserves the Sobol backbone while
borrowing the best part of nonuniform grids: not all tail regions
deserve the same local resolution. The `genzsobol` backend is the first
implemented version that uses reordering and embedded deterministic
stopping together, while `backend = "auto"` adds the outer deterministic
cross-check policy that makes the error attribute more useful in
unattended workflows.

## Accuracy requests and warning

Users can request diagnostic error control with the familiar `abseps`
and `releps` arguments:

``` r
DeterministicQMC(maxpts = 100000, abseps = 1e-6, releps = 0,
                 backend = "auto")
```

    ## $maxpts
    ## [1] 1e+05
    ## 
    ## $abseps
    ## [1] 1e-06
    ## 
    ## $releps
    ## [1] 0
    ## 
    ## $backend
    ## [1] "auto"
    ## 
    ## $start
    ## [1] 0
    ## 
    ## attr(,"class")
    ## [1] "DeterministicQMC" "GenzBretz"

For
[`DeterministicQMC()`](https://keaven.github.io/mvtnorm/reference/algorithms.md),
this is a deterministic stability target. If the reported diagnostic is
larger than the requested absolute or relative tolerance after the
available `maxpts` budget is spent, the result is returned with
`inform = 1` and the message `"Completion with error > abseps"`. This is
useful for unattended deterministic workflows because the caller can
detect that the reproducible diagnostic did not certify the requested
target.

The warning is that this diagnostic is not the same object as the
randomized Monte Carlo error reported by
[`GenzBretz()`](https://keaven.github.io/mvtnorm/reference/algorithms.md).
It is not a confidence interval, and it is not a rigorous deterministic
error bound. It should be read as a reproducible convergence and
stability check.

## Practical interpretation

Use
[`GenzBretz()`](https://keaven.github.io/mvtnorm/reference/algorithms.md)
for general-purpose production calculations where the mature randomized
algorithm and stochastic error assessment are desired. Use
[`Miwa()`](https://keaven.github.io/mvtnorm/reference/algorithms.md) as
a low-dimensional deterministic normal-probability check. Use
[`DeterministicQMC()`](https://keaven.github.io/mvtnorm/reference/algorithms.md)
when reproducibility of the numerical path is more important than
stochastic error assessment. The deterministic `error` attribute should
be read as a diagnostic, not a proof that the true integration error is
above or below the target.
