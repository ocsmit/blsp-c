#ifndef BLSP_H_
#define BLSP_H_
#include "common.h"

typedef struct BLSP_Model {
  gsl_vector *theta_mean_vector;
  gsl_vector *theta_stddev_vector;

  gsl_matrix *theta_hat_matrix;

  gsl_vector *global_log_likelihood_vector;
  gsl_matrix *attempts_matrix;
  gsl_matrix *acceptence_matrix;

  gsl_matrix *metropolis_hastings_matrix;

} BLSP_Model_t;

BLSP_Model_t *BLSP_Model_alloc(size_t year_cnt, size_t obs_cnt,
			       size_t parameter_cnt);


#endif
