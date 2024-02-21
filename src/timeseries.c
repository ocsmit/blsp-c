/**
 *   \file timeseries.c
 *   \brief Methods for creating and working with timeseries
 */

#include "timeseries.h"

TimeSeries_t *TimeSeries_alloc(size_t n_years, size_t n_obs) {
  TimeSeries_t *TSData = malloc(sizeof(TimeSeries_t));
  if (TSData == NULL) {
    printf("Couldn't allocate memory for BLSP_Data\n");
    exit(EXIT_FAILURE);
  }

  TSData->nidx = n_years;
  TSData->nobs = n_obs;

  // right now we are overidding these externally so we just wont allocate them
  TSData->data = gsl_vector_calloc(n_obs);
  TSData->time = gsl_vector_calloc(n_obs);
  TSData->tidx = gsl_vector_calloc(n_years + 1);

  return TSData;
}

TimeSeries_t *TimeSeries_init(size_t n_years, size_t n_obs, gsl_vector *obs,
                              gsl_vector *doy, gsl_vector *year_idx_vec) {
  TimeSeries_t *TSData = malloc(sizeof(TimeSeries_t));
  if (TSData == NULL) {
    printf("Couldn't allocate memory for BLSP_Data\n");
    exit(EXIT_FAILURE);
  }

  TSData->nidx = n_years;
  TSData->nobs = n_obs;

  TSData->data = obs;
  TSData->time = doy;

  return TSData;
}

void TimeSeries_free(TimeSeries_t *TSData) {
  free(TSData->data);
  free(TSData->time);
  free(TSData->tidx);
  free(TSData);
}

unsigned int TimeSeries_year(TimeSeries_t *TSData, size_t i) {
  return gsl_vector_get(TSData->tidx, i);
}

double TimeSeries_doy(TimeSeries_t *TSData, size_t i) {
  return gsl_vector_get(TSData->time, i);
}

gsl_vector_view TimeSeries_doy_year(TimeSeries_t *TSData, size_t year_idx) {
  size_t idx_s, idx_e;

  if (year_idx > (TSData->nidx) - 1)
    exit(EXIT_FAILURE); // TODO: Add message

  idx_s = TimeSeries_year(TSData, year_idx);
  idx_e = TimeSeries_year(TSData, year_idx + 1);

  return gsl_vector_subvector(TSData->time, idx_s, idx_e - idx_s);
}

double TimeSeries_obs(TimeSeries_t *TSData, size_t i) {
  return gsl_vector_get(TSData->data, i);
}

gsl_vector_view TimeSeries_obs_year(TimeSeries_t *TSData, size_t year_idx) {
  size_t idx_s, idx_e;

  if (year_idx > (TSData->nidx) - 1)
    exit(EXIT_FAILURE); // TODO: Add message

  idx_s = TimeSeries_year(TSData, year_idx);
  idx_e = TimeSeries_year(TSData, year_idx + 1);

  return gsl_vector_subvector(TSData->data, idx_s, idx_e - idx_s);
}
