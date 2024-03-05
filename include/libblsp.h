#ifndef LIBBLSP_H_
#define LIBBLSP_H_

#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct BLSP_TimeSeries_T BLSP_TimeSeries_T;
/**
 * @struct blsp_workspace
 *   Holds all required data and structures for BLSP
*/
typedef struct BLSP_Workspace_T BLSP_Workspace_T;

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
extern BLSP_Workspace_T *BLSP_Workspace_alloc(const size_t nobs,
                                                    const size_t nyrs,
                                                    const size_t nburn,
                                                    const size_t niter,
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
extern void BLSP_Workspace_free(BLSP_Workspace_T *w);


/**
 * @brief      Get sample matrix
 *
 * @details    Returns a gsl_matrix* containing all MCMC samples indexed by year and iteration
 *
 * @param      w
 *
 * @return     gsl_matrix *
 */
extern gsl_matrix *BLSP_Workspace_samples(BLSP_Workspace_T *w);


/**
 * @brief      BLSP_TimeSeries_alloc
 *
 * @details    allocates memory for a timeseries struct
 *
 * @param      n_years
 * @param      n_obs
 *
 * @return     *TimeSeries_t
 */
extern BLSP_TimeSeries_T *BLSP_TimeSeries_alloc(size_t n_years, size_t n_obs);

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
extern void BLSP_TimeSeries_free(BLSP_TimeSeries_T *TSData);

extern void BLSP_TimeSeries_set_data(BLSP_TimeSeries_T *TSData, gsl_vector *data);
extern void BLSP_TimeSeries_set_time(BLSP_TimeSeries_T *TSData, gsl_vector *data);
extern void BLSP_TimeSeries_set_tidx(BLSP_TimeSeries_T *TSData, gsl_vector *data);


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
extern int BLSP_sampler(BLSP_TimeSeries_T *X, const gsl_vector *theta_mu,
                 const gsl_vector *theta_sd, BLSP_Workspace_T *w);

#endif // !LIBBLSP_H_
