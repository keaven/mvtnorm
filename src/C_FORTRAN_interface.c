/* $Id: C_FORTRAN_interface.c 357 2020-02-07 14:14:01Z thothorn $
*
*  wrapper for calling R's random number generator from
*  the original FORTRAN code
*
*/

#include "mvtnorm.h"
#include <math.h>
#ifdef __APPLE__
#include <dlfcn.h>
#endif

static int use_dqmc = 0;
static unsigned int dqmc_counter = 0;

void mvtnorm_set_dqmc(int new_use_dqmc, unsigned int start)
{
    use_dqmc = new_use_dqmc;
    dqmc_counter = start;
}

static double dqmc_unifrnd(void)
{
    const double alpha = 0.61803398874989484820; /* 1 / golden ratio */
    double x;

    dqmc_counter++;
    x = fmod(0.5 + alpha * (double) dqmc_counter, 1.0);
    if (x <= 0.0)
        x += alpha;
    if (x >= 1.0)
        x -= floor(x);
    return x;
}

double F77_SUB(unifrnd)(void) { return use_dqmc ? dqmc_unifrnd() : unif_rand(); }
double F77_SUB(sqrtqchisqint)(int *n, double *p) {
    return(sqrt(qchisq(p[0], (double) n[0], 0, 0)));
}
double F77_SUB(phid)(double *x){ return pnorm(*x, 0.0, 1.0, 1, 0); }
double F77_SUB(studnt)(int *nu, double *x){ return pt(x[0], (double) nu[0], 1, 0); }

double F77_SUB(mvphi)(double const *z){
  return pnorm5(*z, 0., 1., 1L, 0L);
}

double F77_SUB(mvphnv)(double const *p){
  return qnorm5(*p, 0., 1., 1L, 0L);
}

SEXP R_dqmc_metal_available(void)
{
    int available = 0;

#ifdef __APPLE__
    void *metal = dlopen("/System/Library/Frameworks/Metal.framework/Metal", RTLD_LAZY);
    if (metal) {
        typedef void *(*create_device_fun)(void);
        create_device_fun create_device =
            (create_device_fun) dlsym(metal, "MTLCreateSystemDefaultDevice");
        if (create_device && create_device())
            available = 1;
        dlclose(metal);
    }
#endif

    return ScalarLogical(available);
}
