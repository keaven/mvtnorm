# Multiple Lower Triangular or Symmetric Matrices

A class representing multiple lower triangular or symmetric matrices and
some methods.

## Usage

``` r
ltMatrices(object, diag = FALSE, byrow = FALSE, names = TRUE)
syMatrices(object, diag = FALSE, byrow = FALSE, names = TRUE)
# S3 method for class 'ltMatrices'
as.array(x, symmetric = FALSE, ...)
# S3 method for class 'syMatrices'
as.array(x, ...)
# S3 method for class 'ltMatrices'
diagonals(x, ...)
# S3 method for class 'syMatrices'
diagonals(x, ...)
# S3 method for class 'matrix'
diagonals(x, ...)
# S3 method for class 'integer'
diagonals(x, ...)
diagonals(x) <- value
# S3 method for class 'ltMatrices'
diagonals(x) <- value
# S3 method for class 'syMatrices'
diagonals(x) <- value
# S3 method for class 'ltMatrices'
solve(a, b, transpose = FALSE, ...)
# S3 method for class 'syMatrices'
chol(x, ...)
# S3 method for class 'chol'
aperm(a, perm, ...)
# S3 method for class 'invchol'
aperm(a, perm, ...)
# S3 method for class 'ltMatrices'
aperm(a, perm, ...)
# S3 method for class 'syMatrices'
aperm(a, perm, ...)
deperma(chol = solve(invchol), permuted_chol = solve(permuted_invchol), 
        invchol, permuted_invchol, perm, score_schol)
# S3 method for class 'ltMatrices'
Mult(x, y, transpose = FALSE, ...)
# S3 method for class 'syMatrices'
Mult(x, y, ...)
Tcrossprod(x, diag_only = FALSE)
Crossprod(x, diag_only = FALSE)
# S3 method for class 'ltMatrices'
tcrossprod(x, y = NULL, ...)
# S3 method for class 'syMatrices'
tcrossprod(x, y = NULL, ...)
# S3 method for class 'ltMatrices'
crossprod(x, y = NULL, ...)
# S3 method for class 'syMatrices'
crossprod(x, y = NULL, ...)
logdet(x)
Lower_tri(x, diag = FALSE, byrow = attr(x, "byrow"))
is.ltMatrices(x)
is.syMatrices(x)
as.ltMatrices(x)
# S3 method for class 'ltMatrices'
as.ltMatrices(x)
# S3 method for class 'syMatrices'
as.ltMatrices(x)
as.syMatrices(x)
is.chol(x)
is.invchol(x)
as.chol(x)
as.invchol(x)
chol2cov(x)
invchol2chol(x)
chol2invchol(x)
invchol2cov(x)
invchol2pre(x)
chol2pre(x)
Dchol(x, D = 1 / sqrt(Tcrossprod(x, diag_only = TRUE)))
invcholD(x, D = sqrt(Tcrossprod(solve(x), diag_only = TRUE)))
chol2cor(x)
invchol2cor(x)
chol2pc(x)
invchol2pc(x)
vectrick(C, S, A, transpose = c(TRUE, TRUE))
standardize(chol, invchol)
destandardize(chol = solve(invchol), invchol, score_schol)
as.ltMatrices(x)
```

## Arguments

- object:

  a `matrix` representing the lower triagular elements of \\N\\ lower
  triangular matrix, each of dimension \\J \times J\\. Dimensions of
  `object` depend on `diag`: With diagonal elements, `object` is a
  \\J(J+1)/2 \times N\\ matrix, otherwise, the number of rows is
  \\J(J - 1) / 2\\.

- diag:

  logical, `object` contains diagonal elements if `TRUE`, otherwise unit
  diagonal elements are assumed.

- byrow:

  logical, `object` represents matrices in row-major order if `TRUE` or,
  otherwise, in column-major order.

- names:

  logical or character vector of length \\J\\.

- symmetric:

  logical, object is interpreted as a symmetric matrix if `TRUE`.

- diag_only:

  logical, compute diagonal elements of crossproduct only if `TRUE`.

- x,chol,invchol,permuted_chol,permuted_invchol:

  object of class `ltMatrices` or `syMatrices` (for `chol`).

- value:

  a matrix of diagonal elements to be assigned (of dimension \\J \times
  N\\).

- a:

  object of class `ltMatrices`.

- perm:

  a permutation of the covariance matrix corresponding to `a`.

- D:

  a matrix (of dimension \\J \times N\\) of diagonal elements to be
  multiplied with.

- y:

  matrix with \\J\\ rows.

- b:

  matrix with \\J\\ rows.

- C:

  an object of class `ltMatrices`.

- S:

  an object of class `ltMatrices` or a matrix with \\J^2\\ rows
  representing multiple \\J x J\\ matrices (columns of vec operators).

- A:

  an object of class `ltMatrices`.

- transpose:

  a logical of length two indicating if `A` or `B` shall be transposed
  in `vectrick`. For `solve`, this argument being true computes
  `solve(t(a), b)` (in absence of a
  [`t()`](https://rdrr.io/r/base/t.html) method for `ltMatrices`
  objects).

- score_schol:

  score matrix for a standardized `chol` object.

- ...:

  additional arguments, currently ignored.

## Details

`ltMatrices` interprets a matrix as lower triangular elements of
multiple lower triangular matrices. The corresponding class can be used
to store such matrices efficiently. Matrix multiplications, solutions to
linear systems, explicite inverses, and crossproducts can be computed
based on such objects. Details can be found in the `lmvnorm_src` package
vignette.

`syMatrices` only store the lower triangular parts of multiple symmetric
matrices.

## Value

The constructor `ltMatrices` returns objects of class `ltMatrices` with
corresponding methods. The constructor `syMatrices` returns objects of
class `syMatrices` with a reduced set of methods.

## See also

[`vignette("lmvnorm_src", package = "mvtnorm")`](https://cran.rstudio.com/web/packages/mvtnorm/vignettes/lmvnorm_src.pdf)

## Examples

``` r
  J <- 4L
  N <- 2L
  dm <- paste0("d", 1:J)
  xm <- paste0("x", 1:N)
  (C <- ltMatrices(matrix(runif(N * J * (J + 1) / 2), 
                          ncol = N, dimnames = list(NULL, xm)), 
                   diag = TRUE, names = dm))
#> , , x1
#> 
#>             d1          d2          d3          d4
#> d1 0.738284444           .           .           .
#> d2 0.129652295 0.143109764           .           .
#> d3 0.465802073 0.510915912 0.078370095           .
#> d4 0.887361601 0.880502964 0.205885270 0.516684529
#> 
#> , , x2
#> 
#>             d1          d2          d3          d4
#> d1 0.005416701           .           .           .
#> d2 0.101740554 0.548209178           .           .
#> d3 0.932764725 0.377695584 0.789203319           .
#> d4 0.161109834 0.979345519 0.824223205 0.211213747
#> 

  ## dimensions and names
  dim(C)
#> [1] 2 4 4
  dimnames(C)
#> [[1]]
#> [1] "x1" "x2"
#> 
#> [[2]]
#> [1] "d1" "d2" "d3" "d4"
#> 
#> [[3]]
#> [1] "d1" "d2" "d3" "d4"
#> 
  names(C)
#>  [1] "d1.d1" "d2.d1" "d3.d1" "d4.d1" "d2.d2" "d3.d2" "d4.d2" "d3.d3" "d4.d3"
#> [10] "d4.d4"

  ## subset
  C[,2:3]
#> , , x1
#> 
#>           d2        d3
#> d2 0.1431098         .
#> d3 0.5109159 0.0783701
#> 
#> , , x2
#> 
#>           d2        d3
#> d2 0.5482092         .
#> d3 0.3776956 0.7892033
#> 

  ## multiplication
  y <- matrix(runif(N * J), nrow = J)
  Mult(C, y)
#>            x1         x2
#> d1 0.09654679 0.00046862
#> d2 0.14911600 0.03084242
#> d3 0.53953880 0.61960378
#> d4 1.35925523 0.74123925
  C 
#> , , x1
#> 
#>             d1          d2          d3          d4
#> d1 0.738284444           .           .           .
#> d2 0.129652295 0.143109764           .           .
#> d3 0.465802073 0.510915912 0.078370095           .
#> d4 0.887361601 0.880502964 0.205885270 0.516684529
#> 
#> , , x2
#> 
#>             d1          d2          d3          d4
#> d1 0.005416701           .           .           .
#> d2 0.101740554 0.548209178           .           .
#> d3 0.932764725 0.377695584 0.789203319           .
#> d4 0.161109834 0.979345519 0.824223205 0.211213747
#> 

  ## solve
  solve(C)
#> , , x1
#> 
#>               d1            d2            d3            d4
#> d1    1.35449149             .             .             .
#> d2   -1.22712054    6.98764342             .             .
#> d3   -0.05065117  -45.55434310   12.75996918             .
#> d4   -0.21485478    6.24428893   -5.08451396    1.93541696
#> 
#> , , x2
#> 
#>               d1            d2            d3            d4
#> d1  184.61420371             .             .             .
#> d2  -34.26201552    1.82412123             .             .
#> d3 -201.79971528   -0.87298484    1.26710060             .
#> d4  805.53088286   -5.05133121   -4.94462949    4.73454032
#> 
  solve(C, y)
#>             x1         x2
#> d1   0.1771293  15.971696
#> d2   6.2925802  -2.890803
#> d3 -40.9691953 -16.652722
#> d4   6.8416089  69.365116

  ## tcrossprod
  Tcrossprod(C)
#> , , x1
#> 
#>            d1         d2        d3        d4
#> d1 0.54506392 0.09572027 0.3438944 0.6551253
#> d2 0.09572027 0.03729012 0.1335094 0.2410570
#> d3 0.34389442 0.13350936 0.4841485 0.8793331
#> d4 0.65512527 0.24105704 0.8793331 1.8720477
#> 
#> , , x2
#> 
#>              d1           d2          d3           d4
#> d1 2.934065e-05 0.0005510982 0.005052508 0.0008726839
#> d2 5.510982e-04 0.3108844429 0.301956185 0.5532776054
#> d3 5.052508e-03 0.3019561851 1.635545864 1.1706517365
#> d4 8.726839e-04 0.5532776054 1.170651736 1.7090291628
#> 
  tcrossprod(C)
#> , , x1
#> 
#>            d1         d2        d3        d4
#> d1 0.54506392 0.09572027 0.3438944 0.6551253
#> d2 0.09572027 0.03729012 0.1335094 0.2410570
#> d3 0.34389442 0.13350936 0.4841485 0.8793331
#> d4 0.65512527 0.24105704 0.8793331 1.8720477
#> 
#> , , x2
#> 
#>              d1           d2          d3           d4
#> d1 2.934065e-05 0.0005510982 0.005052508 0.0008726839
#> d2 5.510982e-04 0.3108844429 0.301956185 0.5532776054
#> d3 5.052508e-03 0.3019561851 1.635545864 1.1706517365
#> d4 8.726839e-04 0.5532776054 1.170651736 1.7090291628
#> 

  ## convert to matrix
  as.array(solve(C[1,]))[,,1]
#>             d1         d2        d3       d4
#> d1  1.35449149   0.000000  0.000000 0.000000
#> d2 -1.22712054   6.987643  0.000000 0.000000
#> d3 -0.05065117 -45.554343 12.759969 0.000000
#> d4 -0.21485478   6.244289 -5.084514 1.935417
```
