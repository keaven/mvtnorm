library("mvtnorm")

corr <- matrix(c(1, .45, .2,
                 .45, 1, -.3,
                 .2, -.3, 1), 3, 3)
lower <- c(-Inf, -.5, -1)
upper <- c(.7, 1.2, .8)
ctrl <- DeterministicQMC(maxpts = 25000, abseps = 0, releps = 0,
                         backend = "c")
sobol <- DeterministicQMC(maxpts = 25000, abseps = 0, releps = 0,
                          backend = "sobol")
hybrid <- DeterministicQMC(maxpts = 25000, abseps = 0, releps = 0,
                           backend = "hybrid")
genzsobol <- DeterministicQMC(maxpts = 25000, abseps = 0, releps = 0,
                              backend = "genzsobol")
auto <- DeterministicQMC(maxpts = 25000, abseps = 1e-5, releps = 0,
                         backend = "auto")
ref <- DeterministicQMC(maxpts = 25000, abseps = 0, releps = 0,
                        backend = "fortran")

set.seed(1)
p1 <- pmvnorm(lower = lower, upper = upper, corr = corr, algorithm = ctrl)
set.seed(99)
p2 <- pmvnorm(lower = lower, upper = upper, corr = corr, algorithm = ctrl)
stopifnot(identical(p1, p2))
stopifnot(abs(p1 - pmvnorm(lower = lower, upper = upper, corr = corr,
                           algorithm = ref)) < 5e-3)

set.seed(1)
ps1 <- pmvnorm(lower = lower, upper = upper, corr = corr, algorithm = sobol)
set.seed(99)
ps2 <- pmvnorm(lower = lower, upper = upper, corr = corr, algorithm = sobol)
stopifnot(identical(ps1, ps2))
stopifnot(abs(ps1 - pmvnorm(lower = lower, upper = upper, corr = corr,
                            algorithm = ref)) < 5e-3)

set.seed(1)
ph1 <- pmvnorm(lower = lower, upper = upper, corr = corr, algorithm = hybrid)
set.seed(99)
ph2 <- pmvnorm(lower = lower, upper = upper, corr = corr, algorithm = hybrid)
stopifnot(identical(ph1, ph2))
stopifnot(abs(ph1 - pmvnorm(lower = lower, upper = upper, corr = corr,
                            algorithm = ref)) < 5e-3)

set.seed(1)
pgs1 <- pmvnorm(lower = lower, upper = upper, corr = corr, algorithm = genzsobol)
set.seed(99)
pgs2 <- pmvnorm(lower = lower, upper = upper, corr = corr, algorithm = genzsobol)
stopifnot(identical(pgs1, pgs2))
stopifnot(abs(pgs1 - pmvnorm(lower = lower, upper = upper, corr = corr,
                             algorithm = ref)) < 5e-3)

set.seed(1)
pa1 <- pmvnorm(lower = lower, upper = upper, corr = corr, algorithm = auto)
set.seed(99)
pa2 <- pmvnorm(lower = lower, upper = upper, corr = corr, algorithm = auto)
stopifnot(identical(pa1, pa2))
stopifnot(identical(attr(attr(pa1, "algorithm"), "backend"),
                    "auto-genzsobol"))
stopifnot(abs(pa1 - pmvnorm(lower = lower, upper = upper, corr = corr,
                            algorithm = ref)) < 5e-3)

set.seed(1)
t1 <- suppressWarnings(pmvt(lower = lower, upper = upper, corr = corr, df = 5,
                            algorithm = ctrl))
set.seed(99)
t2 <- suppressWarnings(pmvt(lower = lower, upper = upper, corr = corr, df = 5,
                            algorithm = ctrl))
stopifnot(identical(t1, t2))

set.seed(1)
ta1 <- pmvt(lower = lower, upper = upper, corr = corr, df = 5,
            algorithm = auto)
set.seed(99)
ta2 <- pmvt(lower = lower, upper = upper, corr = corr, df = 5,
            algorithm = auto)
stopifnot(identical(ta1, ta2))
stopifnot(identical(attr(attr(ta1, "algorithm"), "backend"), "fortran"))

tight <- DeterministicQMC(maxpts = 384, abseps = 1e-12, releps = 0,
                          backend = "auto")
ptight <- pmvnorm(lower = lower, upper = upper, corr = corr,
                  algorithm = tight)
stopifnot(identical(attr(ptight, "msg"),
                    "Completion with error > abseps"))

gb <- GenzBretz(maxpts = 250000, abseps = 0, releps = 0)
set.seed(7)
pgb <- pmvnorm(lower = lower, upper = upper, corr = corr, algorithm = gb)
stopifnot(abs(p1 - pgb) < 5e-4)
stopifnot(abs(ps1 - pgb) < 5e-4)
stopifnot(abs(ph1 - pgb) < 5e-4)
stopifnot(abs(pgs1 - pgb) < 5e-4)
stopifnot(abs(pa1 - pgb) < 5e-4)
