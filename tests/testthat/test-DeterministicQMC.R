dqmc_problem <- function() {
    list(
        corr = matrix(c(1, .45, .2,
                        .45, 1, -.3,
                        .2, -.3, 1), 3, 3),
        lower = c(-Inf, -.5, -1),
        upper = c(.7, 1.2, .8)
    )
}

test_that("deterministic C backends are repeatable and close to Miwa", {
    problem <- dqmc_problem()
    ref <- suppressWarnings(
        pmvnorm(lower = problem$lower, upper = problem$upper,
                corr = problem$corr, algorithm = Miwa(steps = 512))
    )

    for (backend in c("c", "sobol", "hybrid", "genzsobol")) {
        alg <- DeterministicQMC(maxpts = 25000, abseps = 0, releps = 0,
                                backend = backend)
        set.seed(1)
        p1 <- pmvnorm(lower = problem$lower, upper = problem$upper,
                      corr = problem$corr, algorithm = alg)
        set.seed(99)
        p2 <- pmvnorm(lower = problem$lower, upper = problem$upper,
                      corr = problem$corr, algorithm = alg)

        expect_identical(p1, p2)
        expect_lt(abs(as.numeric(p1) - as.numeric(ref)), 1e-3)
        expect_identical(attr(attr(p1, "algorithm"), "backend"), backend)
    }
})

test_that("auto backend chooses deterministic normal path and t fallback", {
    problem <- dqmc_problem()
    alg <- DeterministicQMC(maxpts = 25000, abseps = 1e-5, releps = 0,
                            backend = "auto")

    p <- pmvnorm(lower = problem$lower, upper = problem$upper,
                 corr = problem$corr, algorithm = alg)
    expect_identical(attr(attr(p, "algorithm"), "backend"),
                     "auto-genzsobol")

    tprob <- pmvt(lower = problem$lower, upper = problem$upper,
                  corr = problem$corr, df = 5, algorithm = alg)
    expect_identical(attr(attr(tprob, "algorithm"), "backend"), "fortran")
})

test_that("deterministic tolerance target is reported as a diagnostic", {
    problem <- dqmc_problem()
    alg <- DeterministicQMC(maxpts = 384, abseps = 1e-12, releps = 0,
                            backend = "auto")

    p <- pmvnorm(lower = problem$lower, upper = problem$upper,
                 corr = problem$corr, algorithm = alg)
    expect_identical(attr(p, "msg"), "Completion with error > abseps")
})
