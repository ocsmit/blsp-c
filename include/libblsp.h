#ifndef LIBBLSP_H_
#define LIBBLSP_H_

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdlib.h>

#define PARAMETER_COUNT 7

typedef struct BLSP_TimeSeries_T BLSP_TimeSeries_T;
/**
 * @struct blsp_workspace
 *   Holds all required data and structures for BLSP
 */
typedef struct BLSP_Fit_T BLSP_Fit_T;

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
extern BLSP_Fit_T *BLSP_Fit_alloc(const size_t nobs, const size_t nyrs,
                                  const size_t nburn, const size_t niter,
                                  const size_t npar);


/**
 * @brief      Free workspace memory
 *
 * @details    Frees workspace and all objects needed by MCMC sampler
 *
 * @param      w
 *
 * @return     void
 */
extern void BLSP_Fit_free(BLSP_Fit_T *w);

/**
 * @brief      Get sample matrix
 *
 * @details    Returns a gsl_matrix* containing all MCMC samples indexed by year
 * and iteration
 *
 * @param      w
 *
 * @return     gsl_matrix *
 */
extern gsl_matrix *BLSP_Fit_samples(BLSP_Fit_T *w);


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
BLSP_Fit_T *BLSP_sampler(BLSP_TimeSeries_T *X, const gsl_vector *theta_mu,
                         const gsl_vector *theta_sd, size_t iterations,
                         size_t burn);



BLSP_Fit_T *BLSP_Fit_sample(double *obs_vec, double *doy_vec, size_t obs_count,
                            size_t *idx_vec, size_t idx_count,
                            double init_theta_mu[PARAMETER_COUNT],
                            double init_theta_sd[PARAMETER_COUNT],
                            size_t iterations, size_t burn);

#endif // !LIBBLSP_H_
