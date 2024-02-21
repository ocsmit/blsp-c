#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

SEXP run_blsp(SEXP data_vector, SEXP doy_vector, SEXP year_idx_vector,
              SEXP init_theta_mu, SEXP init_theta_sd, SEXP iterations,
              SEXP burn);
