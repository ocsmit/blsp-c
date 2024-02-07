#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>


#include "../../src/blsp.h"
#include "../../src/blsp_data.h"
#include "../../src/sampler.h"

SEXP generate_data(SEXP data_vector, SEXP doy_vector, SEXP year_idx_vector) {


  int n_obs = length(data_vector);
  int n_years = length(year_idx_vector);


  double *data_vec = REAL(data_vector);
  double *doy_vec = REAL(doy_vector);
  int *year_idx_vec = INTEGER(year_idx_vector);

  SEXP theta_matrix;
  PROTECT(theta_matrix = Rf_allocMatrix(REALSXP, n_years, 7));

  BLSP_Data_t *BLSP_Data = BLSP_Data_alloc(n_years, n_obs);

  gsl_vector_view pp = gsl_vector_view_array(data_vec, n_obs);
  BLSP_Data->obs = &pp.vector;

  gsl_vector_view pv = gsl_vector_view_array(doy_vec, n_obs);
  BLSP_Data->doy = &pv.vector;

  for (size_t i = 1; i < n_years + 1; ++i) {
    gsl_vector_set(BLSP_Data->year_idx, i,
		   gsl_vector_get(BLSP_Data->year_idx, i - 1) +
		       year_idx_vec[i - 1]);
  }


  double theta_mean[7] = {-2.030660, 0.9576184,   119.2568116,
				    8.7140161, 291.0681806, 9.9803594,
				    -6.6876395};

  double theta_sd[7] = {0.03324139, 0.02499047, 0.46365685,
				  0.35526145, 0.75144906, 0.55121990,
				  0.09682220};

  gsl_matrix* theta = blsp_sampler(BLSP_Data, theta_mean, theta_sd, 0.8, 7000, 2000);

  for (int i = 0; i < n_years; i++) {
    for (int p = 0; p < 7; p++) {
      REAL(theta_matrix)[i + n_years*p] = gsl_matrix_get(theta, i, p);
    }
  }

  UNPROTECT(1);
  return theta_matrix;
}
