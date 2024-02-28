#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdio.h>

#include "blsp_r.h"

#include "../../src/blsp.h"
#include "../../src/timeseries.h"
#include "../../src/workspace.h"
#include "R_ext/Print.h"

#define GSL_TO_REALSXP(dst, src, size)                                         \
  memcpy(REAL(dst), src, sizeof(double) * size);

#define GSLVEC_TO_REALSXP(dst, src) GSL_TO_REALSXP(dst, src->data, src->size)
#define GSLMAT_TO_REALSXP(dst, src)                                            \
  GSL_TO_REALSXP(dst, src->data, src->size1 * src->size2)

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

  // NOTE: The following allocates memory for the timeseries to make sure that
  // we have enough memory for it to fit, but I could instead initialize it by
  // assigning the vector views since we know that the data already fits in
  // memory since we passing it from R
  TimeSeries_t *X = TimeSeries_alloc(nyrs, nobs);

  gsl_vector_view data_view = gsl_vector_view_array(data_vec, nobs);
  gsl_vector_memcpy(X->data, &data_view.vector);

  gsl_vector_view doy_view = gsl_vector_view_array(doy_vec, nobs);
  gsl_vector_memcpy(X->time, &doy_view.vector);

  gsl_vector_view idx_view = gsl_vector_view_array(year_idx_vec, nyrs + 1);
  gsl_vector_memcpy(X->tidx, &idx_view.vector);

  // This instead uses pointers to the underlying memory passed from R
  /// TimeSeries_t *Xi = TimeSeries_init(nyrs, nobs, &data_view.vector,
  /// &doy_view.vector, &yr_idx_view.vector);

  blsp_workspace *w = blsp_workspace_sampler_alloc(nobs, nyrs, nburn, niter, 7);

  gsl_vector_view theta_mu_view = gsl_vector_view_array(theta_mu, 7);
  gsl_vector_view theta_sd_view = gsl_vector_view_array(theta_sd, 7);

  int status = blsp_sampler(X, &theta_mu_view.vector, &theta_sd_view.vector, w);

  // Copy over sampling mat
  GSLMAT_TO_REALSXP(theta_matrix, w->parameter_track);

  // Clean up
  TimeSeries_free(X);
  blsp_free_workspace(w);

  // Pop R memory stack
  UNPROTECT(1); // -> theta_matrix

  return theta_matrix;
}
