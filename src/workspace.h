#ifndef BLSP_WORKSPACE_H_
#define BLSP_WORKSPACE_H_

#include "common.h"

/**
 * @struct blsp_workspace
 *   Holds all required data and structures for BLSP
*/
typedef struct {
  size_t nobs;                 /* number of observations */
  size_t nyrs;                 /* number of years */
  size_t npar;                 /* number of parameters */
  gsl_matrix *theta_hat;       /* parameter means for each year */
  gsl_matrix *mh_mat;          /* sampling grid holding current stddev */
  gsl_matrix *accept_mat;      /* tracks candidate acceptences */
  gsl_matrix *attempt_mat;     /* tracks attempts */
  gsl_matrix *parameter_track; /* matrix to track mcmc draws */
  gsl_vector *vi_vec;          /* vegetation index vector */
  gsl_vector *doy_vec;         /* date of year vector */

  // Initialize a vector to hold a /copy/ of the current year's vector so that
  // we can update the years parameter vector with the current parameters draw
  // without overwriting the value in the original theta_hat matrix until we
  // decide to accept the candidate draw.
  gsl_vector *can_p_vec;    /* candidate vector */
  gsl_vector *theta_mu_vec; /* global parameter means */
  gsl_vector *theta_sd_vec; /* global parameter stddevs */
  gsl_vector *loglike_vec;  /* global loglikelihoods */
  double sigma;             /* global variance (noise) */
  size_t niter;             /* number of iterations */
  size_t nburn;             /* number of burnin iterations */
} BLSP_Workspace_T;



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
BLSP_Workspace_T *BLSP_Workspace_alloc(const size_t nobs,
                                             const size_t nyrs,
                                             const size_t nburn,
                                             const size_t niter,
                                             const size_t npar);

/* FIXME: Does not work why? */
/* #define blsp_workspace_alloc(nobs, nyrs, nburn, niter)                         \ */
/*   BLSP_workspace_sampler_alloc(nobs, nyrs, nburn, niter, PARAMETER_COUNT) */

gsl_matrix *BLSP_Workspace_samples(BLSP_Workspace_T *w);

/**
 * @brief      initialize sampling grids
 *
 * @details    initializes sampling grids and starting thetas for MCMC
 *
 * @param      mu
 * @param      sd
 * @param      w
 *
 *
 * @return     void
 */
void BLSP_Workspace_init_thetas(const gsl_vector *mu, const gsl_vector *sd, BLSP_Workspace_T *w);


/**
 * @brief      Free workspace memory
 *
 * @details    Frees workspace and all objects needed by MCMC sampler
 *
 * @param      w
 *
 * @return     void
 */
void BLSP_Workspace_free(BLSP_Workspace_T *w);

/**
 * @brief      Fill gsl_matrix row wise
 *
 * @details    Fills all values in gsl_matrix with a vector row wise
 *
 * @param      m
 * @param      x
 *
 * @return     void
 */
void fill_matrix_by_row(gsl_matrix *m, gsl_vector *x);

#endif //! BLSP_WORKSPACE_H_
