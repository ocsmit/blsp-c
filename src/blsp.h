#ifndef SAMPLER_H_
#define SAMPLER_H_

#include "blsp_math.h"
#include "timeseries.h"
#include "workspace.h"

#define sample_matrix gsl_matrix

typedef struct {
  // BLSP_TimeSeries_T *data;
  sample_matrix *samples;
  struct series {
    gsl_vector *data;
    gsl_vector *time;
    gsl_vector *tidx;
  };
  size_t len; // => nidx
  size_t cnt; // => nobs
  // gsl_matrix *quantiles; TODO
} BLSP_T;

// MCMC_T;

/**
 * @brief      BLSP Algorithm
 *
 * @details    Runs MCMC sampler for BLSP
 *
 * @param      X          TimeSeries_t
 * @param      theta_mu   mean parameter values to initialize the parameters for
 * each year
 * @param      theta_sd   standard deviation parameter values to initialize the
 * parameters
 * @param      w          workspace
 *
 * @return     int
 */
/* int BLSP_sampler(BLSP_TimeSeries_T *X, const gsl_vector *theta_mu, */
/*                  const gsl_vector *theta_sd, BLSP_Workspace_T *w); */

BLSP_Fit_T *BLSP_sampler(BLSP_TimeSeries_T *X, const gsl_vector *theta_mu,
                         const gsl_vector *theta_sd, size_t iterations,
                         size_t burn);

int BLSP_mcmc(double *obs_vec, double *doy_vec, size_t obs_count, size_t *idx_vec,
              size_t idx_count, double init_theta_mu[7],
              double init_theta_sd[7], size_t iterations, size_t burn);
// void BLSP_sample(BLSP_Fit_T *X);

#endif
