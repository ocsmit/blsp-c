#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

#define GSL_TO_REALSXP(dst, src, size)                                         \
  memcpy(REAL(dst), src, sizeof(double) * size);

#define GSLVEC_TO_REALSXP(dst, src) GSL_TO_REALSXP(dst, src->data, src->size)

#define GSLMAT_TO_REALSXP(dst, src)                                            \
  GSL_TO_REALSXP(dst, src->data, src->size1 * src->size2)

/**
 * @brief      Run BLSP algorithm
 *
 * @details    Handles R data type conversions for C blsp algorithm.
 *
 * @param      data_vector          vector of observation data
 * @param      doy_vector           vector of DOY's which correspond to each
 *                                  observation
 * @param      year_idx_vector      vector of which contains the number of
 *                                  observations contained within each year
 * @param      init_theta_mu        vector of initial parameter means
 * @param      init_theta_sd        vector of initial parameter std devs
 * @param      iterations           number of iterations to sample
 * @param      burn                 number of iterations for burn-in
 *
 * @return     SEXP
 */
SEXP run_blsp(SEXP data_vector, SEXP doy_vector, SEXP year_idx_vector,
              SEXP init_theta_mu, SEXP init_theta_sd, SEXP iterations,
              SEXP burn);
