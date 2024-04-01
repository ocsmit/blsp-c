/**
 *   \file timeseries.c
 *   \brief Methods for creating and working with timeseries
 */

#include "timeseries.h"

#define view_to_vec(v) &v.vector

BLSP_TimeSeries_T *BLSP_TimeSeries_alloc(size_t n_years, size_t n_obs) {
  BLSP_TimeSeries_T *TSData = malloc(sizeof(BLSP_TimeSeries_T));
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

void BLSP_TimeSeries_init(BLSP_TimeSeries_T *X, double *obs, double *doy,
                                        double *year_idx_vec) {

  gsl_vector_view obs_view, doy_view, idx_view;

  obs_view = gsl_vector_view_array(obs, X->nobs);
  doy_view = gsl_vector_view_array(doy, X->nobs);
  idx_view = gsl_vector_view_array(year_idx_vec, X->nidx + 1);

  gsl_vector_memcpy(X->data, view_to_vec(obs_view));
  gsl_vector_memcpy(X->time, view_to_vec(doy_view));
  gsl_vector_memcpy(X->tidx, view_to_vec(idx_view));
}

void BLSP_TimeSeries_free(BLSP_TimeSeries_T *TSData) {
  free(TSData->data);
  free(TSData->time);
  free(TSData->tidx);
  free(TSData);
}

void BLSP_TimeSeries_set_data(BLSP_TimeSeries_T *TSData, gsl_vector *data) {
  gsl_vector_memcpy(TSData->data, data);
}
void BLSP_TimeSeries_set_time(BLSP_TimeSeries_T *TSData, gsl_vector *data) {
  gsl_vector_memcpy(TSData->time, data);
}
void BLSP_TimeSeries_set_tidx(BLSP_TimeSeries_T *TSData, gsl_vector *data) {
  gsl_vector_memcpy(TSData->tidx, data);
}

unsigned int BLSP_TimeSeries_year(BLSP_TimeSeries_T *TSData, size_t i) {
  return gsl_vector_get(TSData->tidx, i);
}

double BLSP_TimeSeries_doy(BLSP_TimeSeries_T *TSData, size_t i) {
  return gsl_vector_get(TSData->time, i);
}

gsl_vector_view BLSP_TimeSeries_doy_year(BLSP_TimeSeries_T *TSData,
                                         size_t year_idx) {
  size_t idx_s, idx_e;

  if (year_idx > (TSData->nidx) - 1)
    exit(EXIT_FAILURE); // TODO: Add message

  idx_s = BLSP_TimeSeries_year(TSData, year_idx);
  idx_e = BLSP_TimeSeries_year(TSData, year_idx + 1);

  return gsl_vector_subvector(TSData->time, idx_s, idx_e - idx_s);
}

double BLSP_TimeSeries_obs(BLSP_TimeSeries_T *TSData, size_t i) {
  return gsl_vector_get(TSData->data, i);
}

gsl_vector_view BLSP_TimeSeries_obs_year(BLSP_TimeSeries_T *TSData,
                                         size_t year_idx) {
  size_t idx_s, idx_e;

  if (year_idx > (TSData->nidx) - 1)
    exit(EXIT_FAILURE); // TODO: Add message

  idx_s = BLSP_TimeSeries_year(TSData, year_idx);
  idx_e = BLSP_TimeSeries_year(TSData, year_idx + 1);

  return gsl_vector_subvector(TSData->data, idx_s, idx_e - idx_s);
}
