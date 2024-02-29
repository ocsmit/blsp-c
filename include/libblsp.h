#ifndef LIBBLSP_H_
#define LIBBLSP_H_

#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct TimeSeries_t TimeSeries_t;
typedef struct blsp_workspace blsp_workspace;

/**
 * @brief      Allocate workspace
 *
 * @details    Allocates memory for workspace
 *
 * @param      nobs
 * @param      nyrs
 * @param      nburn
 * @param      niter
 * @param      npar
 *
 * @return     blsp_workspace
 */
extern blsp_workspace *blsp_workspace_sampler_alloc(const size_t nobs,
                                                    const size_t nyrs,
                                                    const size_t nburn,
                                                    const size_t niter,
                                                    const size_t npar);

extern void blsp_free_workspace(blsp_workspace *w);

extern gsl_matrix *blsp_workspace_samples(blsp_workspace *w);

extern TimeSeries_t *TimeSeries_alloc(size_t n_years, size_t n_obs);
extern void TimeSeries_free(TimeSeries_t *TSData);
extern void TimeSeries_set_data(TimeSeries_t *TSData, gsl_vector *data);
extern void TimeSeries_set_time(TimeSeries_t *TSData, gsl_vector *data);
extern void TimeSeries_set_tidx(TimeSeries_t *TSData, gsl_vector *data);


extern int blsp_sampler(TimeSeries_t *X, const gsl_vector *theta_mu,
                 const gsl_vector *theta_sd, blsp_workspace *w);

#endif // !LIBBLSP_H_
