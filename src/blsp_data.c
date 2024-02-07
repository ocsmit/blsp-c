#include "blsp_data.h"

BLSP_Data_t *BLSP_Data_alloc(size_t n_years, size_t n_obs) {
  BLSP_Data_t *BLSP_Data = malloc(sizeof(BLSP_Data_t));
  if (BLSP_Data == NULL) {
    printf("Couldn't allocate memory for BLSP_Data\n");
    exit(EXIT_FAILURE);
  }

  BLSP_Data->n_years = n_years;
  BLSP_Data->n_obs = n_obs;

  // right now we are overidding these externally so we just wont allocate them
  // BLSP_Data->obs = gsl_vector_calloc(n_obs);
  // BLSP_Data->doy = gsl_vector_calloc(n_obs);
  BLSP_Data->year_idx = gsl_vector_calloc(n_years + 1);

  return BLSP_Data;
}


void BLSP_Data_free(BLSP_Data_t* BLSP_Data) {
  free(BLSP_Data->obs);
  free(BLSP_Data->doy);
  free(BLSP_Data->year_idx);
  free(BLSP_Data);
}

unsigned int BLSP_Data_get_year(BLSP_Data_t *BLSP_Data, size_t i) {
  return gsl_vector_get(BLSP_Data->year_idx, i);
}


double BLSP_Data_get_doy(BLSP_Data_t *BLSP_Data, size_t i) {
  return gsl_vector_get(BLSP_Data->doy, i);
}

gsl_vector_view BLSP_Data_get_year_doy(BLSP_Data_t *BLSP_Data, size_t year_idx) {
  size_t idx_s, idx_e;

  if (year_idx > (BLSP_Data->n_years) - 1)
    exit(EXIT_FAILURE); // TODO: Add message

  idx_s = BLSP_Data_get_year(BLSP_Data, year_idx);
  idx_e = BLSP_Data_get_year(BLSP_Data, year_idx + 1);

  return gsl_vector_subvector(BLSP_Data->doy, idx_s, idx_e - idx_s);
}

double BLSP_Data_get_obs(BLSP_Data_t *BLSP_Data, size_t i) {
  return gsl_vector_get(BLSP_Data->obs, i);
}

gsl_vector_view BLSP_Data_get_year_obs(BLSP_Data_t *BLSP_Data, size_t year_idx) {
  size_t idx_s, idx_e;

  if (year_idx > (BLSP_Data->n_years) - 1)
    exit(EXIT_FAILURE); // TODO: Add message

  idx_s = BLSP_Data_get_year(BLSP_Data, year_idx);
  idx_e = BLSP_Data_get_year(BLSP_Data, year_idx + 1);

  return gsl_vector_subvector(BLSP_Data->obs, idx_s, idx_e - idx_s);
}
