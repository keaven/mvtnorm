#include "mvtnorm.h"
#include <float.h>
#include <math.h>

#define DQMC_HALTON 0
#define DQMC_SOBOL 1
#define DQMC_HYBRID 2
#define DQMC_GENZSOBOL 3
#define SOBOL_BITS 32
#define SOBOL_MAX_DIM 32
#define HYBRID_MIX 0.10
#define HYBRID_POWER 2.0
#define GENZSOBOL_MIN_CHECK 64
#define TWO32 4294967296.0

typedef struct {
    int degree;
    unsigned int coefficient;
    unsigned int m[7];
} sobol_param;

static const sobol_param sobol_params[SOBOL_MAX_DIM] = {
    {0, 0, {0, 0, 0, 0, 0, 0, 0}},
    {1, 0, {1, 0, 0, 0, 0, 0, 0}},
    {2, 1, {1, 3, 0, 0, 0, 0, 0}},
    {3, 1, {1, 3, 1, 0, 0, 0, 0}},
    {3, 2, {1, 1, 1, 0, 0, 0, 0}},
    {4, 1, {1, 1, 3, 3, 0, 0, 0}},
    {4, 4, {1, 3, 5, 13, 0, 0, 0}},
    {5, 2, {1, 1, 5, 5, 17, 0, 0}},
    {5, 4, {1, 1, 5, 5, 5, 0, 0}},
    {5, 7, {1, 1, 7, 11, 19, 0, 0}},
    {5, 11, {1, 1, 5, 1, 1, 0, 0}},
    {5, 13, {1, 1, 1, 3, 11, 0, 0}},
    {5, 14, {1, 3, 5, 5, 31, 0, 0}},
    {6, 1, {1, 3, 3, 9, 7, 49, 0}},
    {6, 13, {1, 1, 1, 15, 21, 21, 0}},
    {6, 16, {1, 3, 1, 13, 27, 49, 0}},
    {6, 19, {1, 1, 1, 15, 7, 5, 0}},
    {6, 22, {1, 3, 1, 15, 13, 25, 0}},
    {6, 25, {1, 1, 5, 5, 19, 61, 0}},
    {7, 1, {1, 3, 7, 11, 23, 15, 103}},
    {7, 4, {1, 3, 7, 13, 13, 15, 69}},
    {7, 7, {1, 1, 3, 13, 7, 35, 63}},
    {7, 8, {1, 3, 5, 9, 1, 25, 53}},
    {7, 14, {1, 3, 1, 13, 9, 35, 107}},
    {7, 19, {1, 3, 1, 5, 27, 61, 31}},
    {7, 21, {1, 1, 5, 11, 19, 41, 61}},
    {7, 28, {1, 3, 5, 3, 3, 13, 69}},
    {7, 31, {1, 1, 7, 13, 1, 19, 1}},
    {7, 32, {1, 3, 7, 5, 13, 19, 59}},
    {7, 37, {1, 1, 3, 9, 25, 29, 41}},
    {7, 41, {1, 3, 5, 13, 23, 1, 55}},
    {7, 42, {1, 3, 7, 3, 13, 59, 17}}
};

static double corr_at(double *corr, int n, int i, int j)
{
    return corr[i + j * n];
}

static int chol_lower(double *corr, double *L, int n)
{
    double sum, lij, tol = sqrt(DBL_EPSILON);

    for (int i = 0; i < n * n; i++)
        L[i] = 0.0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            sum = corr_at(corr, n, i, j);
            for (int k = 0; k < j; k++)
                sum -= L[i * n + k] * L[j * n + k];
            if (i == j) {
                if (sum <= tol)
                    return 3;
                L[i * n + j] = sqrt(sum);
            } else {
                lij = L[j * n + j];
                if (fabs(lij) <= tol)
                    return 3;
                L[i * n + j] = sum / lij;
            }
        }
    }
    return 0;
}

static int next_prime(int x)
{
    int is_prime;

    for (;;) {
        x++;
        is_prime = x >= 2;
        for (int d = 2; d * d <= x && is_prime; d++)
            if (x % d == 0)
                is_prime = 0;
        if (is_prime)
            return x;
    }
}

static double halton(unsigned int index, int base)
{
    double f = 1.0 / (double) base;
    double r = 0.0;

    while (index > 0) {
        r += f * (double) (index % (unsigned int) base);
        index /= (unsigned int) base;
        f /= (double) base;
    }
    return r;
}

static int sobol_init_dir(unsigned int *directions, int ndim)
{
    if (ndim > SOBOL_MAX_DIM)
        return 4;

    for (int dim = 0; dim < ndim; dim++) {
        unsigned int *v = directions + (size_t) dim * SOBOL_BITS;
        const sobol_param *param = sobol_params + dim;
        int degree = param->degree;
        unsigned int coefficient = param->coefficient;

        if (dim == 0) {
            for (int bit = 0; bit < SOBOL_BITS; bit++)
                v[bit] = 1u << (SOBOL_BITS - 1 - bit);
            continue;
        }

        for (int bit = 0; bit < degree; bit++)
            v[bit] = param->m[bit] << (SOBOL_BITS - 1 - bit);

        for (int bit = degree; bit < SOBOL_BITS; bit++) {
            unsigned int next = v[bit - degree] ^ (v[bit - degree] >> degree);
            for (int k = 1; k < degree; k++)
                if ((coefficient >> (degree - 1 - k)) & 1u)
                    next ^= v[bit - k];
            v[bit] = next;
        }
    }
    return 0;
}

static unsigned int sobol_uint(unsigned int index, int coord,
                               const unsigned int *directions)
{
    unsigned int gray = index ^ (index >> 1);
    unsigned int x = 0u;
    const unsigned int *v = directions + (size_t) coord * SOBOL_BITS;

    for (int bit = 0; gray && bit < SOBOL_BITS; bit++, gray >>= 1)
        if (gray & 1u)
            x ^= v[bit];

    return x;
}

static double sobol_unit(unsigned int x)
{
    return ((double) x + 0.5) / TWO32;
}

static double tail_warp(double u, double *jacobian)
{
    double s = 2.0 * u - 1.0;
    double a = fabs(s);
    double tail, dtail;

    if (a == 0.0) {
        tail = 0.5;
        dtail = 0.0;
    } else {
        tail = 0.5 + 0.5 * (s > 0.0 ? 1.0 : -1.0) * pow(a, HYBRID_POWER);
        dtail = HYBRID_POWER * pow(a, HYBRID_POWER - 1.0);
    }

    *jacobian = (1.0 - HYBRID_MIX) + HYBRID_MIX * dtail;
    return (1.0 - HYBRID_MIX) * u + HYBRID_MIX * tail;
}

static int sobol_update_bit(unsigned int index)
{
    int bit = 0;

    while (index && !(index & 1u) && bit < SOBOL_BITS) {
        index >>= 1;
        bit++;
    }
    return bit;
}

static double clamp_prob(double x)
{
    if (x < DBL_EPSILON)
        return DBL_EPSILON;
    if (x > 1.0 - DBL_EPSILON)
        return 1.0 - DBL_EPSILON;
    return x;
}

static void prob_bounds(int infin, double lower, double upper,
                        double shift, double scale, double *d, double *e)
{
    if (infin < 0) {
        *d = 0.0;
        *e = 1.0;
    } else if (infin == 0) {
        *d = 0.0;
        *e = pnorm((upper - shift) / scale, 0.0, 1.0, 1, 0);
    } else if (infin == 1) {
        *d = pnorm((lower - shift) / scale, 0.0, 1.0, 1, 0);
        *e = 1.0;
    } else {
        *d = pnorm((lower - shift) / scale, 0.0, 1.0, 1, 0);
        *e = pnorm((upper - shift) / scale, 0.0, 1.0, 1, 0);
    }
}

static void init_identity_order(int *order, int n)
{
    for (int i = 0; i < n; i++)
        order[i] = i;
}

static double marginal_width(int infin, double lower, double upper,
                             double delta, double variance)
{
    double d, e;

    prob_bounds(infin, lower - delta, upper - delta, 0.0, sqrt(variance),
                &d, &e);
    return e - d;
}

static void genzsobol_order(double *lower, double *upper, double *corr,
                            double *delta, int *infin, int *order, int n)
{
    double *width = (double *) R_alloc((size_t) n, sizeof(double));

    init_identity_order(order, n);
    for (int i = 0; i < n; i++)
        width[i] = marginal_width(infin[i], lower[i], upper[i], delta[i],
                                  corr_at(corr, n, i, i));

    for (int i = 0; i < n - 1; i++) {
        int best = i;
        for (int j = i + 1; j < n; j++)
            if (width[order[j]] < width[order[best]])
                best = j;
        if (best != i) {
            int tmp = order[i];
            order[i] = order[best];
            order[best] = tmp;
        }
    }
}

static void reorder_problem(double *lower, double *upper, double *corr,
                            double *delta, int *infin, int *order, int n,
                            double **rlower, double **rupper, double **rcorr,
                            double **rdelta, int **rinfin)
{
    double *lower2 = (double *) R_alloc((size_t) n, sizeof(double));
    double *upper2 = (double *) R_alloc((size_t) n, sizeof(double));
    double *corr2 = (double *) R_alloc((size_t) n * (size_t) n,
                                       sizeof(double));
    double *delta2 = (double *) R_alloc((size_t) n, sizeof(double));
    int *infin2 = (int *) R_alloc((size_t) n, sizeof(int));

    for (int i = 0; i < n; i++) {
        int oi = order[i];
        lower2[i] = lower[oi];
        upper2[i] = upper[oi];
        delta2[i] = delta[oi];
        infin2[i] = infin[oi];
        for (int j = 0; j < n; j++)
            corr2[i + j * n] = corr[oi + order[j] * n];
    }

    *rlower = lower2;
    *rupper = upper2;
    *rcorr = corr2;
    *rdelta = delta2;
    *rinfin = infin2;
}

static SEXP dqmc_mvn(SEXP lower, SEXP upper, SEXP corr, SEXP delta, SEXP infin,
                     SEXP maxpts, SEXP abseps, SEXP releps, SEXP start,
                     int rule)
{
    int n = LENGTH(lower);
    int M = INTEGER(maxpts)[0];
    int start_value = INTEGER(start)[0];
    int inform = 0;
    int nblocks = 8;
    int adaptive = (rule == DQMC_GENZSOBOL);
    double value = 0.0, error = 0.0;
    double eps_abs = (abseps == R_NilValue) ? -1.0 : REAL(abseps)[0];
    double eps_rel = (releps == R_NilValue) ? 0.0 : REAL(releps)[0];
    SEXP ans, names;

    if (M < 1)
        M = 1;
    if (nblocks > M)
        nblocks = M;

    PROTECT(ans = allocVector(VECSXP, 3));
    PROTECT(names = allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, mkChar("value"));
    SET_STRING_ELT(names, 1, mkChar("error"));
    SET_STRING_ELT(names, 2, mkChar("inform"));
    setAttrib(ans, R_NamesSymbol, names);

    double *dlower = REAL(lower);
    double *dupper = REAL(upper);
    double *dcorr = REAL(corr);
    double *ddelta = REAL(delta);
    int *iinfin = INTEGER(infin);

    if (rule == DQMC_GENZSOBOL && n > 1) {
        int *order = (int *) R_alloc((size_t) n, sizeof(int));
        genzsobol_order(dlower, dupper, dcorr, ddelta, iinfin, order, n);
        reorder_problem(dlower, dupper, dcorr, ddelta, iinfin, order, n,
                        &dlower, &dupper, &dcorr, &ddelta, &iinfin);
    }

    if (n == 1) {
        double d, e;
        prob_bounds(iinfin[0], dlower[0] - ddelta[0], dupper[0] - ddelta[0],
                    0.0, sqrt(corr_at(dcorr, n, 0, 0)), &d, &e);
        value = e - d;
    } else {
        double *L = (double *) R_alloc((size_t) n * (size_t) n, sizeof(double));
        double *z = (double *) R_alloc((size_t) (n - 1), sizeof(double));
        double *block_sum = (double *) R_alloc((size_t) nblocks, sizeof(double));
        int *block_n = (int *) R_alloc((size_t) nblocks, sizeof(int));
        int *bases = NULL;
        unsigned int *directions = NULL;
        unsigned int *sobol_state = NULL;

        inform = chol_lower(dcorr, L, n);
        if (inform == 0) {
            if (rule == DQMC_HALTON) {
                int p = 1;
                bases = (int *) R_alloc((size_t) (n - 1), sizeof(int));
                for (int j = 0; j < n - 1; j++) {
                    p = next_prime(p);
                    bases[j] = p;
                }
            } else {
                directions = (unsigned int *) R_alloc((size_t) (n - 1) *
                                                      SOBOL_BITS,
                                                      sizeof(unsigned int));
                inform = sobol_init_dir(directions, n - 1);
                if (inform == 0)
                    sobol_state = (unsigned int *) R_alloc((size_t) (n - 1),
                                                           sizeof(unsigned int));
            }
            if (inform == 0) {
                unsigned int idx =
                    (unsigned int) ((start_value < 0 ? 0 : start_value) + 1);

                for (int b = 0; b < nblocks; b++) {
                    block_sum[b] = 0.0;
                    block_n[b] = 0;
                }

                if (rule != DQMC_HALTON)
                    for (int j = 0; j < n - 1; j++)
                        sobol_state[j] = sobol_uint(idx, j, directions);

                int count = 0, stop = 0, has_prev = 0;
                int next_check = GENZSOBOL_MIN_CHECK;
                double prev_mean = 0.0, stability_error = 0.0;

                for (int m = 0; m < M && !stop; m++) {
                    double d, e, emd, f, q, shift, u, jacobian;
                    int block = (int) (((double) m * (double) nblocks) / (double) M);
                    if (block >= nblocks)
                        block = nblocks - 1;

                    prob_bounds(iinfin[0], dlower[0] - ddelta[0],
                                dupper[0] - ddelta[0], 0.0, L[0], &d, &e);
                    emd = e - d;
                    f = emd;

                    for (int j = 1; j < n; j++) {
                        u = (rule != DQMC_HALTON) ?
                            sobol_unit(sobol_state[j - 1]) :
                            halton(idx, bases[j - 1]);
                        jacobian = 1.0;
                        if (rule == DQMC_HYBRID)
                            u = tail_warp(u, &jacobian);
                        q = d + u * emd;
                        z[j - 1] = qnorm(clamp_prob(q), 0.0, 1.0, 1, 0);
                        shift = ddelta[j];
                        for (int k = 0; k < j; k++)
                            shift += L[j * n + k] * z[k];
                        prob_bounds(iinfin[j], dlower[j], dupper[j],
                                    shift, L[j * n + j], &d, &e);
                        emd = e - d;
                        f *= emd * jacobian;
                    }
                    value += f;
                    count++;
                    block_sum[block] += f;
                    block_n[block]++;

                    if (adaptive && (count == next_check || count == M)) {
                        double mean = value / (double) count;
                        if (has_prev) {
                            double tol = eps_abs < 0.0 ? 0.0 : eps_abs;
                            double rtol = eps_rel * fabs(mean);
                            if (rtol > tol)
                                tol = rtol;
                            stability_error = fabs(mean - prev_mean);
                            if (count >= GENZSOBOL_MIN_CHECK &&
                                stability_error <= tol)
                                stop = 1;
                        }
                        prev_mean = mean;
                        has_prev = 1;
                        if (next_check < M)
                            next_check = next_check > M / 2 ? M : 2 * next_check;
                    }

                    if (!stop && m + 1 < M) {
                        idx++;
                        if (rule != DQMC_HALTON) {
                            int bit = sobol_update_bit(idx);
                            if (bit < SOBOL_BITS) {
                                for (int j = 0; j < n - 1; j++)
                                    sobol_state[j] ^=
                                        directions[(size_t) j * SOBOL_BITS + bit];
                            } else {
                                for (int j = 0; j < n - 1; j++)
                                    sobol_state[j] = sobol_uint(idx, j, directions);
                            }
                        }
                    }
                }
                value /= (double) count;

                if (adaptive) {
                    error = stability_error;
                } else if (nblocks > 1) {
                    double mean = value, var = 0.0;
                    for (int b = 0; b < nblocks; b++) {
                        double bm = block_sum[b] / (double) block_n[b];
                        var += (bm - mean) * (bm - mean);
                    }
                    var /= (double) (nblocks - 1);
                    error = 3.0 * sqrt(var / (double) nblocks);
                }
            }
        }
    }

    SET_VECTOR_ELT(ans, 0, ScalarReal(value));
    SET_VECTOR_ELT(ans, 1, ScalarReal(error));
    SET_VECTOR_ELT(ans, 2, ScalarInteger(inform));
    UNPROTECT(2);
    return ans;
}

SEXP R_dqmc_mvn(SEXP lower, SEXP upper, SEXP corr, SEXP delta, SEXP infin,
                SEXP maxpts, SEXP start)
{
    return dqmc_mvn(lower, upper, corr, delta, infin, maxpts,
                    R_NilValue, R_NilValue, start, DQMC_HALTON);
}

SEXP R_dqmc_mvn_sobol(SEXP lower, SEXP upper, SEXP corr, SEXP delta,
                      SEXP infin, SEXP maxpts, SEXP start)
{
    return dqmc_mvn(lower, upper, corr, delta, infin, maxpts,
                    R_NilValue, R_NilValue, start, DQMC_SOBOL);
}

SEXP R_dqmc_mvn_hybrid(SEXP lower, SEXP upper, SEXP corr, SEXP delta,
                       SEXP infin, SEXP maxpts, SEXP start)
{
    return dqmc_mvn(lower, upper, corr, delta, infin, maxpts,
                    R_NilValue, R_NilValue, start, DQMC_HYBRID);
}

SEXP R_dqmc_mvn_genzsobol(SEXP lower, SEXP upper, SEXP corr, SEXP delta,
                          SEXP infin, SEXP maxpts, SEXP abseps,
                          SEXP releps, SEXP start)
{
    return dqmc_mvn(lower, upper, corr, delta, infin, maxpts,
                    abseps, releps, start, DQMC_GENZSOBOL);
}
