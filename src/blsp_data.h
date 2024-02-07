#ifndef BLSP_DATA_H_
#define BLSP_DATA_H_

#include "common.h"

/*
 * Input data format
 * - For now two vector blocks - 1 for signal data and 1 for doy vector
 *    - future: Block of memory := sizeof(double) * obs_count * 2
 *      - First chunk of memory will be the signal data and second chunk
 *        will be DOY vector
 * - year_idx, this will give index that points to the first value of that
 *   year. The size will be n_year + 1
 *
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


// Observation
typedef struct BLSP_Data {
  unsigned int n_years; // Total number of years
  unsigned int n_obs;   // total number of observations
  gsl_vector *obs;      // GSL vector holding signal data
  gsl_vector *doy;      // GLS vector holding DOY for each observation
  gsl_vector *year_idx; // GSL vector holding idx for each year
} BLSP_Data_t;

BLSP_Data_t *BLSP_Data_alloc(size_t n_years, size_t n_obs);
void BLSP_Data_free(BLSP_Data_t* BLSP_Data);

unsigned int BLSP_Data_get_year(BLSP_Data_t *BLSP_Data, size_t i);

double BLSP_Data_get_doy(BLSP_Data_t *BLSP_Data, size_t i);
gsl_vector_view BLSP_Data_get_year_doy(BLSP_Data_t *BLSP_Data, size_t year_idx);


double BLSP_Data_get_obs(BLSP_Data_t *BLSP_Data, size_t i);
gsl_vector_view BLSP_Data_get_year_obs(BLSP_Data_t *BLSP_Data, size_t year_idx);


#endif // !TS_DATA_H_
