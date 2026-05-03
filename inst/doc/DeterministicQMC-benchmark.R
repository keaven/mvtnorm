### R code from vignette source '/Users/keaven/GitHub/mvtnorm/vignettes/DeterministicQMC-benchmark.Rnw'

###################################################
### code chunk number 1: setup
###################################################
library("mvtnorm")


###################################################
### code chunk number 2: benchmark
###################################################
make_problem <- function(n) {
    list(corr = toeplitz(0.35^(0:(n - 1))),
         lower = rep(-Inf, n),
         upper = seq(-0.2, 0.6, length.out = n))
}

methods <- list(
    list(name = "Default GenzBretz",
         algorithm = GenzBretz(),
         seed = 11L),
    list(name = "Miwa steps=128",
         algorithm = Miwa(steps = 128),
         seed = NULL),
    list(name = "Original randomized QMC",
         algorithm = GenzBretz(maxpts = 25000, abseps = 0, releps = 0),
         seed = 11L),
    list(name = "DeterministicQMC C",
         algorithm = DeterministicQMC(maxpts = 25000, abseps = 0,
                                      releps = 0, backend = "c"),
         seed = NULL),
    list(name = "DeterministicQMC Sobol",
         algorithm = DeterministicQMC(maxpts = 25000, abseps = 0,
                                      releps = 0, backend = "sobol"),
         seed = NULL),
    list(name = "DeterministicQMC Hybrid",
         algorithm = DeterministicQMC(maxpts = 25000, abseps = 0,
                                      releps = 0, backend = "hybrid"),
         seed = NULL),
    list(name = "DeterministicQMC GenzSobol",
         algorithm = DeterministicQMC(maxpts = 25000, abseps = 1e-5,
                                      releps = 0, backend = "genzsobol"),
         seed = NULL),
    list(name = "DeterministicQMC Auto",
         algorithm = DeterministicQMC(maxpts = 50000, abseps = 1e-5,
                                      releps = 0, backend = "auto"),
         seed = NULL)
)

timed_value <- function(problem, method, reps = 3L) {
    times <- numeric(reps)
    value <- NA_real_
    diagnostic <- NA_real_
    target_requested <- is.list(method$algorithm) &&
        ((length(method$algorithm$abseps) == 1L &&
          !is.na(method$algorithm$abseps) && method$algorithm$abseps > 0) ||
         (length(method$algorithm$releps) == 1L &&
          !is.na(method$algorithm$releps) && method$algorithm$releps > 0))
    target_status <- "not requested"
    for (i in seq_len(reps)) {
        if (!is.null(method$seed))
            set.seed(method$seed)
        invisible(gc(FALSE))
        tm <- system.time({
            ans <- pmvnorm(lower = problem$lower, upper = problem$upper,
                           corr = problem$corr,
                           algorithm = method$algorithm)
        })
        times[i] <- unname(tm["elapsed"])
        if (i == 1L) {
            value <- as.numeric(ans)
            diagnostic <- attr(ans, "error")
            if (target_requested)
                target_status <- if (identical(attr(ans, "msg"),
                                               "Normal Completion"))
                    "met" else "not met"
        }
    }
    list(value = value, diagnostic = diagnostic,
         target_status = target_status, elapsed = median(times))
}

benchmark_dimension <- function(n) {
    problem <- make_problem(n)
    reference <- as.numeric(pmvnorm(lower = problem$lower,
                                    upper = problem$upper,
                                    corr = problem$corr,
                                    algorithm = Miwa(steps = 512)))
    rows <- lapply(methods, function(method) {
        ans <- timed_value(problem, method)
        data.frame(dimension = n,
                   method = method$name,
                   elapsed_sec = sprintf("%.4f", ans$elapsed),
                   value = sprintf("%.8g", ans$value),
                   diagnostic = sprintf("%.3e", ans$diagnostic),
                   target_status = ans$target_status,
                   rel_error = sprintf("%.5g",
                       if (reference == 0) NA_real_ else
                           abs(ans$value - reference) / abs(reference)),
                   stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

results <- do.call(rbind, lapply(c(4L, 8L, 12L), benchmark_dimension))
print(results, row.names = FALSE)


