#include "ts_data.h"
#include <stdio.h>
#include <stdlib.h>

TS_Data_t *TS_Data_alloc(size_t n_years, size_t n_obs) {
  TS_Data_t *TS_Data = malloc(sizeof(TS_Data_t));
  if (TS_Data == NULL) {
    printf("Couldn't allocate memory for TS_Data\n");
    exit(EXIT_FAILURE);
  }

  TS_Data->n_years = n_years;
  TS_Data->n_obs = n_obs;


  // right now we are overidding these externally so we just wont allocate them
  // TS_Data->obs = gsl_vector_calloc(n_obs);
  // TS_Data->doy = gsl_vector_calloc(n_obs);
  TS_Data->year_idx = gsl_vector_calloc(n_years + 1);

  return TS_Data;
}


void TS_Data_free(TS_Data_t* TS_Data) {
  free(TS_Data->obs);
  free(TS_Data->doy);
  free(TS_Data->year_idx);
  free(TS_Data);
}

unsigned int TS_Data_get_year(TS_Data_t *TS_Data, size_t i) {
  return gsl_vector_get(TS_Data->year_idx, i);
}


double TS_Data_get_doy(TS_Data_t *TS_Data, size_t i) {
  return gsl_vector_get(TS_Data->doy, i);
}

gsl_vector_view TS_Data_get_year_doy(TS_Data_t *TS_Data, size_t year_idx) {
  size_t idx_s, idx_e;

  if (year_idx > (TS_Data->n_years) - 1)
    exit(EXIT_FAILURE); // TODO: Add message

  idx_s = TS_Data_get_year(TS_Data, year_idx);
  idx_e = TS_Data_get_year(TS_Data, year_idx + 1);

  return gsl_vector_subvector(TS_Data->doy, idx_s, idx_e - idx_s);
}

double TS_Data_get_obs(TS_Data_t *TS_Data, size_t i) {
  return gsl_vector_get(TS_Data->obs, i);
}

gsl_vector_view TS_Data_get_year_obs(TS_Data_t *TS_Data, size_t year_idx) {
  size_t idx_s, idx_e;

  if (year_idx > (TS_Data->n_years) - 1)
    exit(EXIT_FAILURE); // TODO: Add message

  idx_s = TS_Data_get_year(TS_Data, year_idx);
  idx_e = TS_Data_get_year(TS_Data, year_idx + 1);

  return gsl_vector_subvector(TS_Data->obs, idx_s, idx_e - idx_s);
}

