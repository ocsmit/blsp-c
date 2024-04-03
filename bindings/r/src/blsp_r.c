#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <libblsp.h>

#include "blsp_r.h"

#include <R_ext/Print.h>
#include <gsl/gsl_vector.h>

#define view_to_vec(v) &v.vector

SEXP run_blsp(SEXP data_vector, SEXP doy_vector, SEXP year_idx_vector,
              SEXP init_theta_mu, SEXP init_theta_sd, SEXP iterations,
              SEXP burn) {

  int nobs = length(data_vector);
  int nyrs = length(year_idx_vector) - 1; // yr idx vector is + 1
  int niter = asInteger(iterations);
  int nburn = asInteger(burn);

  double *data_vec = REAL(data_vector);
  double *doy_vec = REAL(doy_vector);
  double *year_idx_vec = REAL(year_idx_vector);
  double *theta_mu = REAL(init_theta_mu);
  double *theta_sd = REAL(init_theta_sd);

  SEXP theta_matrix;
  PROTECT(theta_matrix = Rf_allocVector(REALSXP, nyrs * niter * 9));

  BLSP_Fit_T *blsp_fit =
      BLSP_Fit_sample(data_vec, doy_vec, nobs, (size_t *)year_idx_vec, nyrs,
                      theta_mu, theta_sd, niter, nburn);

  // Copy over sampling mat
  GSLMAT_TO_REALSXP(theta_matrix, BLSP_Fit_samples(blsp_fit));

  // Clean up
  // BLSP_TimeSeries_free(X);
  BLSP_Fit_free(blsp_fit);

  // Pop R memory stack
  UNPROTECT(1); // -> theta_matrix

  return theta_matrix;
}
