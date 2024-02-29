// #define DEBUG
#include "blsp.h"
#include "blsp_math.h"
#include "common.h"
#include "timeseries.h"

#define check_file(fptr)                                                       \
  if (!fptr) {                                                                 \
    printf("Error! Could not open file\n");                                    \
    exit(-1);                                                                  \
  }

void fprint_gslmat(const gsl_matrix *m) {
  for (size_t i = 0; i < m->size1; ++i) {
    for (size_t j = 0; j < m->size2; ++j)
      fprintf(stderr, "%f ", gsl_matrix_get(m, i, j));
    fprintf(stderr, "\n");
  }
}

void print_gslmat(const gsl_matrix *m) {
  for (size_t i = 0; i < m->size1; ++i) {
    for (size_t j = 0; j < m->size2; ++j)
      printf("%f ", gsl_matrix_get(m, i, j));
    printf("\n");
  }
}

void fprint_gslvec(const gsl_vector *m) {
  for (size_t i = 0; i < m->size; ++i) {
    fprintf(stderr, "%f ", gsl_vector_get(m, i));
  }
  fprintf(stderr, "\n");
}

// MAIN ***********************************************************************
int main(int argc, char const *argv[]) {

  const int N_PARAMETER = 7;

  // Parameters taken from simulated data in R
  double theta_mean[7] = {-1.998392,  0.960355, 120.702350, 9.263498,
                          288.853856, 9.166685, -6.592421};

  double theta_sd[7] = {0.07057906, 0.05609551, 1.08944966, 0.88183154,
                        1.55979462, 1.20727157, 0.19881890};

  // Copy over theta mean and theta standard deviation to GSL_Vectors
  // Could also just create a GSL array view to the two arrays but would rather
  // just copy the data for now
  gsl_vector *theta_sd_view = gsl_vector_calloc(N_PARAMETER);
  gsl_vector *theta_mean_view = gsl_vector_calloc(N_PARAMETER);
  for (size_t i = 0; i < N_PARAMETER; ++i) {
    gsl_vector_set(theta_mean_view, i, theta_mean[i]);
    gsl_vector_set(theta_sd_view, i, theta_sd[i]);
  }

  // We simulated 20 years and 334 obs
  size_t num_years = 20, num_obs = 334, curr_yr_id = 0;
  double prior_strength = 0.1; // used when updating MCMC

  // GSL vectors to store data
  gsl_vector *doy_vector = gsl_vector_alloc(num_obs);
  gsl_vector *y_vector = gsl_vector_alloc(num_obs);
  gsl_vector *yr_id_vector = gsl_vector_alloc(num_years);

  // Read in test data. GLS will throw error if file not found
  //FILE *file_doy = fopen("./data/doy_full.dat", "r");
  FILE *file_doy = fopen(argv[1], "r");
  check_file(file_doy);
  gsl_vector_fscanf(file_doy, doy_vector);
  fclose(file_doy);

  //FILE *file_y = fopen("./data/y_full.dat", "r");
  FILE *file_y = fopen(argv[2], "r");
  check_file(file_y);
  gsl_vector_fscanf(file_y, y_vector);
  fclose(file_y);

  //FILE *file_yr_id = fopen("./data/year_idx.dat", "r");
  FILE *file_yr_id = fopen(argv[3], "r");
  check_file(file_yr_id);
  gsl_vector_fscanf(file_y, yr_id_vector);
  fclose(file_yr_id);

  TimeSeries_t *X = TimeSeries_alloc(num_years, num_obs);

  gsl_vector_memcpy(X->data, y_vector);
  gsl_vector_memcpy(X->time, doy_vector);

  /* gsl_vector_view idx_view = gsl_vector_view_array(year_idx_vec, nyrs + 1);
   */
  /* gsl_vector_memcpy(X->tidx, &idx_view.vector); */

  // Read in test data. GLS will throw error if file not found
  for (size_t i = 1; i < X->nidx + 1; ++i) {
    gsl_vector_set(X->tidx, i,
                   gsl_vector_get(X->tidx, i - 1) +
                       gsl_vector_get(yr_id_vector, i - 1));
  }

  // This instead uses pointers to the underlying memory passed from R
  /// TimeSeries_t *Xi = TimeSeries_init(nyrs, nobs, &data_view.vector,
  /// &doy_view.vector, &yr_idx_view.vector);

  blsp_workspace *w =
      blsp_workspace_sampler_alloc(num_obs, num_years, 2000, 7000, 7);

  blsp_sampler(X, theta_mean_view, theta_sd_view, w);

  TimeSeries_free(X);
  gsl_vector_free(theta_mean_view);
  gsl_vector_free(theta_sd_view);
  printf("Done!\n");
  return 0;
}
