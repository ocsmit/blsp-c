#ifndef SAMPLER_H_
#define SAMPLER_H_

#include "blsp_math.h"
#include "timeseries.h"
#include "workspace.h"


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
int BLSP_sampler(BLSP_TimeSeries_T *X, const gsl_vector *theta_mu,
                 const gsl_vector *theta_sd, BLSP_Workspace_T *w);


#endif
