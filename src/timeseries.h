/**
 *   \file timeseries.h
 *   \brief data structures and methods for timeseries
 */
#ifndef TIMESERIES_H_
#define TIMESERIES_H_

#include "common.h"

#define TimeSeries_slice gsl_vector_view

/**
 * @struct TimeSeries_t
 * Struct providing interface to timeseries data indexed by year
 *
 * The timeseries struct contains 3 vectors, vals which holds the data, time
 * which holds the temporal info, and idx which provides the start index for how
 * they are grouped.
 *
 *   ┌───┬───┬───┬───┬───┐
 *   │   │   │   │   │   │ year vector
 *   └───┴───┴───┴───┴───┘
 *     │   │   │   │   │
 *   ┌─┘  ┌┘  ┌┘  ┌┘  ┌┘
 *   ▼    ▼   ▼   ▼   ▼
 *   0   12  34  52  79─────────┐ each value is the starting idx of data for
 *   │    │   │   │             │ that year
 *  ┌┘    │   └─┐ └───────┐     └────────┐
 *  │     │     │         │              │
 *  ▼─────▼─────▼─────────▼──────────────▼─────────┐
 *  │                                              │ signal | obs vector
 *  └──────────────────────────────────────────────┘
 *
 *

 */
typedef struct TimeSeries {
  size_t nidx;      /* Total number of years */
  size_t nobs;      /* total number of observations */
  gsl_vector *data; /* GSL vector holding signal data */
  gsl_vector *time; /* GLS vector holding time for each observation */
  gsl_vector *tidx; /* GSL vector holding idx for time*/
} BLSP_TimeSeries;

/* TODO: Implement view/slice mechanism for timeseries */
typedef struct TimeSeries_Slice {
  gsl_vector_view *data;
  gsl_vector_view *time;
} BLSP_TimeSeries_Slice;

/**
 * @brief      TimeSeries_alloc
 *
 * @details    allocates memory for a timeseries struct
 *
 * @param      n_years
 * @param      n_obs
 *
 * @return     *TimeSeries_t
 */
BLSP_TimeSeries *BLSP_TimeSeries_alloc(size_t n_years, size_t n_obs);

/**
 * @brief      Frees TimeSeries_t memory
 *
 * @details    Frees memory allocated by TimeSeries_alloc for a TimeSeries_t
 * object
 *
 * @param      TSData
 *
 * @return     void
 */
void BLSP_TimeSeries_free(BLSP_TimeSeries *TSData);

/* FIXME: segfaults when used with R interface */
BLSP_TimeSeries *BLSP_TimeSeries_init(size_t n_years, size_t n_obs, gsl_vector *obs,
                              gsl_vector *doy, gsl_vector *year_idx_vec);

/**
 * @brief      Get tidx value
 *
 * @details    returns ith location in tidx
 *
 * @param      TSData
 * @param      i
 *
 * @return     unsigned int
 */
unsigned int BLSP_TimeSeries_year(BLSP_TimeSeries *TSData, size_t i);

/**
 * @brief      get time value
 *
 * @details    returns the time of observation i
 *
 * @param      TSData
 * @param      i
 *
 * @return     return type
 */
double BLSP_TimeSeries_doy(BLSP_TimeSeries *TSData, size_t i);

/**
 * @brief      Subset by index
 *
 * @details    Subsets time by the grouping idx
 *
 * @param      TSData
 * @param      i
 *
 * @return     return type
 */
TimeSeries_slice BLSP_TimeSeries_doy_year(BLSP_TimeSeries *TSData, size_t tidx);

/**
 * @brief      get obs value
 *
 * @details    returns the value of observation i
 *
 * @param      TSData
 * @param      i
 *
 * @return     return type
 */
double BLSP_TimeSeries_obs(BLSP_TimeSeries *TSData, size_t i);


/**
 * @brief      Subset by index
 *
 * @details    Subsets time by the grouping idx
 *
 * @param      TSData
 * @param      i
 *
 * @return     return type
 */
TimeSeries_slice BLSP_TimeSeries_obs_year(BLSP_TimeSeries *BLSP_Data, size_t year_idx);

void BLSP_TimeSeries_set_data(BLSP_TimeSeries *TSData, gsl_vector *data);
void BLSP_TimeSeries_set_time(BLSP_TimeSeries *TSData, gsl_vector *data);
void BLSP_TimeSeries_set_tidx(BLSP_TimeSeries *TSData, gsl_vector *data);

#endif // !TIMESERIES_H_
