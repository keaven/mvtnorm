### R code from vignette source '/Users/keaven/GitHub/mvtnorm/vignettes/DeterministicQMC-methods.Rnw'

###################################################
### code chunk number 1: setup
###################################################
library("mvtnorm")


###################################################
### code chunk number 2: constructors
###################################################
GenzBretz()
DeterministicQMC(backend = "auto")
DeterministicQMC(backend = "c")
DeterministicQMC(backend = "sobol")
DeterministicQMC(backend = "hybrid")
DeterministicQMC(backend = "genzsobol")
DeterministicQMC(backend = "fortran")
Miwa()


###################################################
### code chunk number 3: accuracy-request
###################################################
DeterministicQMC(maxpts = 100000, abseps = 1e-6, releps = 0,
                 backend = "auto")


