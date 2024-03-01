#ifndef LIBBLSP_H_
#define LIBBLSP_H_

#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct BLSP_TimeSeries BLSP_TimeSeries;
typedef struct BLSP_Workspace BLSP_Workspace;

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
extern BLSP_Workspace *BLSP_Workspace_alloc(const size_t nobs,
                                                    const size_t nyrs,
                                                    const size_t nburn,
                                                    const size_t niter,
                                                    const size_t npar);

extern void BLSP_Workspace_free(BLSP_Workspace *w);

extern gsl_matrix *BLSP_Workspace_samples(BLSP_Workspace *w);

extern BLSP_TimeSeries *BLSP_TimeSeries_alloc(size_t n_years, size_t n_obs);
extern void BLSP_TimeSeries_free(BLSP_TimeSeries *TSData);
extern void BLSP_TimeSeries_set_data(BLSP_TimeSeries *TSData, gsl_vector *data);
extern void BLSP_TimeSeries_set_time(BLSP_TimeSeries *TSData, gsl_vector *data);
extern void BLSP_TimeSeries_set_tidx(BLSP_TimeSeries *TSData, gsl_vector *data);


extern int BLSP_sampler(BLSP_TimeSeries *X, const gsl_vector *theta_mu,
                 const gsl_vector *theta_sd, BLSP_Workspace *w);

#endif // !LIBBLSP_H_
