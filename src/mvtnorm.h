
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>

void F77_NAME(mvtdst)(int *n, int *nu, double *lower, double *upper,
                      int *infin, double *corr, double *delta,
                      int *maxpts, double *abseps, double *releps,
                      double *error, double *value, int *inform);

void F77_NAME(tvtlrcall)(int *NU, double *H, double *R, double *EPSI, double *TVTL);
void C_tvtlr            (int *NU, double *H, double *R, double *EPSI, double *TVTL);

void F77_NAME(bvtlrcall)(int *NU, double *DH, double *DK, double *R, double *BVTL);
void C_bvtlr            (int *NU, double *DH, double *DK, double *R, double *BVTL);

void C_mvtdst(int *n, int *nu, double *lower, double *upper,
              int *infin, double *corr, double *delta,
              int *maxpts, double *abseps, double *releps,
              double *error, double *value, int *inform, int *rnd);
void C_mvtdst_dqmc(int *n, int *nu, double *lower, double *upper,
                   int *infin, double *corr, double *delta,
                   int *maxpts, double *abseps, double *releps,
                   double *error, double *value, int *inform, int *start);

void mvtnorm_set_dqmc(int use_dqmc, unsigned int start);

extern SEXP R_miwa(SEXP steps, SEXP corr, SEXP upper, SEXP lower, SEXP infin);
extern SEXP R_dqmc_metal_available(void);
extern SEXP R_dqmc_mvn(SEXP lower, SEXP upper, SEXP corr, SEXP delta, SEXP infin,
                       SEXP maxpts, SEXP start);
extern SEXP R_dqmc_mvn_sobol(SEXP lower, SEXP upper, SEXP corr, SEXP delta,
                             SEXP infin, SEXP maxpts, SEXP start);
extern SEXP R_dqmc_mvn_hybrid(SEXP lower, SEXP upper, SEXP corr, SEXP delta,
                              SEXP infin, SEXP maxpts, SEXP start);
extern SEXP R_dqmc_mvn_genzsobol(SEXP lower, SEXP upper, SEXP corr,
                                 SEXP delta, SEXP infin, SEXP maxpts,
                                 SEXP abseps, SEXP releps, SEXP start);
extern SEXP R_ltMatrices_solve (SEXP C, SEXP y, SEXP N, SEXP J, SEXP diag, SEXP transpose);
extern SEXP R_ltMatrices_solve_C (SEXP C, SEXP N, SEXP J, SEXP diag, SEXP transpose);
extern SEXP R_ltMatrices_tcrossprod (SEXP C, SEXP N, SEXP J, SEXP diag, SEXP diag_only, SEXP transpose);
extern SEXP R_ltMatrices_logdet (SEXP C, SEXP N, SEXP J, SEXP diag, SEXP byrow);
extern SEXP R_ltMatrices_Mult (SEXP C, SEXP y, SEXP N, SEXP J, SEXP diag);
extern SEXP R_ltMatrices_Mult_transpose (SEXP C, SEXP y, SEXP N, SEXP J, SEXP diag);
extern SEXP R_ltMatrices_colSumsdnorm (SEXP z, SEXP N, SEXP J);
extern SEXP R_lpmvnorm(SEXP a, SEXP b, SEXP C, SEXP center, SEXP N, SEXP J, SEXP W, SEXP M, SEXP tol, SEXP logLik, SEXP fast);
extern SEXP R_slpmvnorm(SEXP a, SEXP b, SEXP C, SEXP center, SEXP N, SEXP J, SEXP W, SEXP M, SEXP tol, SEXP fast);
extern SEXP R_vectrick(SEXP C, SEXP N, SEXP J, SEXP S, SEXP D, SEXP diag, SEXP trans);
extern SEXP R_syMatrices_chol(SEXP Sigma, SEXP N, SEXP J);
