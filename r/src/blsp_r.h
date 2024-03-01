#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

#define GSL_TO_REALSXP(dst, src, size)                                         \
  memcpy(REAL(dst), src, sizeof(double) * size);

#define GSLVEC_TO_REALSXP(dst, src) GSL_TO_REALSXP(dst, src->data, src->size)

#define GSLMAT_TO_REALSXP(dst, src)                                            \
  GSL_TO_REALSXP(dst, src->data, src->size1 * src->size2)

SEXP run_blsp(SEXP data_vector, SEXP doy_vector, SEXP year_idx_vector,
              SEXP init_theta_mu, SEXP init_theta_sd, SEXP iterations,
              SEXP burn);
