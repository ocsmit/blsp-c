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

  BLSP_TimeSeries_T *X = BLSP_TimeSeries_alloc(nyrs, nobs);
  BLSP_TimeSeries_init(X, data_vec, doy_vec, year_idx_vec);

  gsl_vector_view theta_mu_view = gsl_vector_view_array(theta_mu, 7);
  gsl_vector_view theta_sd_view = gsl_vector_view_array(theta_sd, 7);

  BLSP_Fit_T *fit = BLSP_sampler(X, view_to_vec(theta_mu_view),
                                 view_to_vec(theta_sd_view), niter, nburn);

  // Copy over sampling mat
  GSLMAT_TO_REALSXP(theta_matrix, BLSP_Fit_samples(fit));

  // Clean up
  BLSP_TimeSeries_free(X);
  BLSP_Fit_free(fit);

  // Pop R memory stack
  UNPROTECT(1); // -> theta_matrix

  return theta_matrix;
}
