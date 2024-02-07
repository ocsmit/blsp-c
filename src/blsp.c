#include "blsp.h"

BLSP_Model_t *BLSP_Model_alloc(size_t year_cnt, size_t obs_cnt,
			       size_t parameter_cnt) {

  BLSP_Model_t *BLSP_Model = malloc(sizeof(BLSP_Model));

  if (BLSP_Model == NULL) {
    printf("Couldn't allocate memory for BLSP_Model\n");
    exit(EXIT_FAILURE);
  }

  BLSP_Model->theta_mean_vector = gsl_vector_calloc(parameter_cnt);
  BLSP_Model->theta_stddev_vector = gsl_vector_calloc(parameter_cnt);
  BLSP_Model->global_log_likelihood_vector = gsl_vector_calloc(parameter_cnt);

  BLSP_Model->theta_hat_matrix = gsl_matrix_calloc(year_cnt, parameter_cnt);
  BLSP_Model->attempts_matrix = gsl_matrix_calloc(year_cnt, parameter_cnt);
  BLSP_Model->acceptence_matrix = gsl_matrix_calloc(year_cnt, parameter_cnt);
  BLSP_Model->metropolis_hastings_matrix =
      gsl_matrix_calloc(year_cnt, parameter_cnt);

  // NOTE: Surely there is a way to allocate all the needed memory in a single
  // gsl_block for the vectors and matricies at once?
  /* size_t memory_block_size = */
  /*     ((year_cnt * parameter_cnt) * 4) + (parameter_cnt * 3); */

  /* gsl_block *memory_block = gsl_block_calloc(memory_block_size); */

  /* size_t stride = sizeof(double); */

  /* gsl_vector_alloc_from_block(memory_block, 0, parameter_cnt, stride); */

  return BLSP_Model;
}
