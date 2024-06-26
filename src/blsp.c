#include "blsp.h"
#include "common.h"
#include "timeseries.h"
#include "workspace.h"

#define matrix_incr(X, i, j, y)                                                \
  gsl_matrix_set(X, i, j, gsl_matrix_get(X, i, j) + y);

#define matrix_incr_1(X, i, j)                                                 \
  gsl_matrix_set(X, i, j, gsl_matrix_get(X, i, j) + 1.0);

#define matrix_reset_val(X, i, j) gsl_matrix_set(X, i, j, 0.1);

gsl_vector *vector_from_array(unsigned int n, double *x) {
  gsl_vector *v = gsl_vector_alloc(n);
  unsigned int i;
  for (i = 0; i < n; i++)
    gsl_vector_set(v, i, x[i]);
  return v;
}

BLSP_Fit_T *BLSP_sampler(BLSP_TimeSeries_T *X, const gsl_vector *theta_mu,
                         const gsl_vector *theta_sd, size_t iterations,
                         size_t burn) { //, BLSP_Workspace_T *w) {

  BLSP_Fit_T *w = BLSP_Fit_alloc(X->nobs, X->nidx, burn, iterations, 7);

  // Initialize thetas in workspace
  BLSP_Fit_init_thetas(theta_mu, theta_sd, w);

  double theta_prior_m = 0.1; /* strength of priors */
  double theta_prior_w =
      w->nyrs / (w->nyrs + theta_prior_m); /* coeff for post Create */

  // norm random number generator MOVE THIS TO WORKSPACE STRUCT
  const gsl_rng_type *RNG_t;
  gsl_rng *RNG_ptr;
  gsl_rng_env_setup();
  RNG_t = gsl_rng_taus;
  RNG_ptr = gsl_rng_alloc(RNG_t);
  gsl_rng_set(RNG_ptr, 1);

  // INITIAL LOG LIKELIHOOD
  // Initialize empty (zeroed) vector to hold log likelihoods
  for (size_t year = 0; year < w->nyrs; ++year) {
    // Variable to store summed log likelihood for each year
    double year_i_log_likelihood = 0;

    // Get value vectors for year i

    gsl_vector_view year_i_obs = BLSP_TimeSeries_obs_year(X, year);
    gsl_vector_view year_i_doy = BLSP_TimeSeries_doy_year(X, year);

    // Get parameter vector for year i
    gsl_vector_view curr_theta_hat = gsl_matrix_row(w->theta_hat, year);

    double ll = vi_log_likelihood(&year_i_obs.vector, &year_i_doy.vector,
                                  &curr_theta_hat.vector,
                                  year_i_obs.vector.size, w->sigma);

    // Set log likelihood
    gsl_vector_set(w->loglike_vec, year, ll);
  }

  // Set up metropolis hastings grid and initialize
  // This stores the current standard deviation for sampling a specific
  // parameter in a specific year.
  /* for (size_t j = 0; j < w->nyrs; ++j) { */
  /*   gsl_matrix_set_row(w->mh_mat, j, w->theta_sd_vec); */
  /* } */
  gsl_matrix_scale(w->mh_mat, 5.0);

  /* theta_hat_candidate_vector = */
  /*     gsl_vector_calloc(PARAMETER_COUNT); // Candidate vector */

  for (size_t iter = 0; iter < w->niter; ++iter) {
    for (size_t year = 0; year < w->nyrs; ++year) {

      // Create vector views for the data in the current year
      gsl_vector_view year_i_obs = BLSP_TimeSeries_obs_year(X, year);
      gsl_vector_view year_i_doy = BLSP_TimeSeries_doy_year(X, year);
      gsl_vector_view theta_hat_year = gsl_matrix_row(w->theta_hat, year);

      for (size_t param = 0; param < PARAMETER_COUNT; ++param) {

        // Increment attempt matrix
        matrix_incr_1(w->attempt_mat, year, param);

        // memcopy so we can change the current parameter without setting it in
        // the original
        gsl_vector_memcpy(w->can_p_vec, &theta_hat_year.vector);

        // Draw a candidate for parameter from gaussian distribution
        double candidate_draw = draw_ran_gaussian(
            RNG_ptr, gsl_matrix_get(w->theta_hat, year, param),
            gsl_matrix_get(w->mh_mat, year, param));

        // Update the candidate parameter vector
        gsl_vector_set(w->can_p_vec, param, candidate_draw);

        // Compute the candidate log likelihood for current year
        double candidate_year_i_log_likelihood =
            vi_log_likelihood(&year_i_obs.vector, &year_i_doy.vector,
                              w->can_p_vec, year_i_doy.vector.size, w->sigma);

        // Compute acceptance ratio
        double R = log_likelihood_ratio(
            candidate_draw, candidate_year_i_log_likelihood,
            gsl_matrix_get(w->theta_hat, year, param),
            gsl_vector_get(w->loglike_vec, year),
            gsl_vector_get(w->theta_mu_vec, param),
            gsl_vector_get(w->theta_sd_vec, param));

        // Decide to accept or reject
        double U = logl(gsl_ran_flat(RNG_ptr, 0, 1));

        if (U < R) {
          gsl_matrix_set(w->theta_hat, year, param, candidate_draw);
          gsl_vector_set(w->loglike_vec, year, candidate_year_i_log_likelihood);
          matrix_incr_1(w->accept_mat, year, param);
        }

        // Tune metropolis hastings samples
        if (iter < w->nburn) {
          double attempt_cnt = gsl_matrix_get(w->attempt_mat, year, param);
          double acceptence_cnt = gsl_matrix_get(w->accept_mat, year, param);

          // Tune every 50 iterations
          if (attempt_cnt > BURN_STEPS) {
            double current_sample_stddev =
                gsl_matrix_get(w->mh_mat, year, param);

            // If acceptence rate is low, tighten stddev
            if (acceptence_cnt / attempt_cnt < 0.2) {
              gsl_matrix_set(w->mh_mat, year, param,
                             0.8 * current_sample_stddev);
            }

            // If acceptence rate is high, expand sttdev
            if (acceptence_cnt / attempt_cnt > 0.6) {
              gsl_matrix_set(w->mh_mat, year, param,
                             1.2 * current_sample_stddev);
            }

            // Reset the attempt and acceptence value for current parameter
            matrix_reset_val(w->accept_mat, year, param);
            matrix_reset_val(w->attempt_mat, year, param);
          }
        }
      } // End metropolis hastings for parameter vector

      gsl_vector_view current_th = gsl_matrix_row(w->theta_hat, year);

      // TODO: This could be more efficient
      size_t idx = iter + (year * w->niter);
      gsl_matrix_set(w->parameter_track, idx, 1, iter);
      gsl_matrix_set(w->parameter_track, idx, 0, year);
      for (size_t i = 0; i < 7; ++i) {
        gsl_matrix_set(w->parameter_track, idx, i + 2,
                       gsl_vector_get(&current_th.vector, i));
      }

    } // End year loop

    // **** GIBBS *************************************************************
    // Sample for global mean and stddev

    // THETA MEAN *************************************************************
    // Update the vector storing global mean for each parameter
    //
    // NOTE: When refactoring, see if it's faster to use BLAS to calculate
    // column means instead of iterating fully over each column:
    // gsl_blas_dgemv(CblasTrans, 1.0 / theta_hat->size1, theta_hat,
    // scale_vector, 0.0, col_means);

    for (size_t p = 0; p < PARAMETER_COUNT; ++p) {
      double column_mean = 0, column_sum = 0;
      double theta_mean_draw = 0;

      // calculate column mean for parameter p
      for (size_t year = 0; year < w->nyrs; ++year) {
        column_sum += gsl_matrix_get(w->theta_hat, year, p);
      }
      column_mean = column_sum / w->nyrs;

      theta_mean_draw = draw_ran_gaussian(RNG_ptr,
                                          theta_prior_w * column_mean +
                                              (1 - theta_prior_w) *
                                                  gsl_vector_get(theta_mu, p),
                                          gsl_vector_get(w->theta_sd_vec, p) /
                                              sqrtl((w->nyrs + theta_prior_m)));

      gsl_vector_set(w->theta_mu_vec, p, theta_mean_draw);
    } // END THETA MEAN UPDATE

    // THETA STDDEV ***********************************************************
    // Update the vector storing global stddev for each parameter
    // 1. Calculate SSE
    // 2. Draw new stddev from gamma dist

    for (size_t p = 0; p < PARAMETER_COUNT; ++p) {
      // 1. SSE for parameter p
      double xi = 0, xj = 0, stddev_sse = 0;
      double sse_i = 0;
      double sse_2 = 0;
      for (size_t year = 0; year < w->nyrs; ++year) {
        xi = gsl_matrix_get(w->theta_hat, year, p);
        xj = gsl_vector_get(w->theta_mu_vec, p); // move this out (future)
        sse_i = xi - xj;
        sse_2 = gsl_pow_2(sse_i);

        stddev_sse += sse_2;
      }

      // 2. Draw new stddev from gamma dist
      double gamma_a = 0, gamma_b = 0, theta_stddev_p = 0;
      // Compute parameters
      gamma_a = w->nyrs / 2.0 + 0.1;
      gamma_b = stddev_sse / 2.0 + 0.1;

      // Draw sample
      theta_stddev_p = gsl_ran_gamma(RNG_ptr, gamma_a, 1.0 / gamma_b);

      // intermediate variable for debugging
      double theta_stddev_p_final = 1 / sqrtl(theta_stddev_p);

      gsl_vector_set(w->theta_sd_vec, p, theta_stddev_p_final);
    } // END THETA STDDEV UPDATE

    // GLOBAL STDDEV **********************************************************
    // update overall standard dev (noise)
    double sse_sim = 0;
    for (size_t year = 0; year < w->nyrs; ++year) {
      // 1. Calculate VI SSE for each year and add to SSE

      // Get value vectors for year i
      gsl_vector_view year_i_obs = BLSP_TimeSeries_obs_year(X, year);
      gsl_vector_view year_i_doy = BLSP_TimeSeries_doy_year(X, year);

      // We use the scale vector as a temporary place to memcpy our parameter
      // vector for year i
      // gsl_matrix_get_row(candidate_vector, theta_hat, year);

      gsl_vector_view curr_theta_hat = gsl_matrix_row(w->theta_hat, year);
      // Iterate over each obs date in year i summing up SSE with each i
      double vi_sim = 0;
      for (size_t obs_i = 0; obs_i < year_i_obs.vector.size; ++obs_i) {
        // Calculate the VI for obs_i with year_i parameter vector
        double doy_i = gsl_vector_get(&year_i_doy.vector, obs_i);
        vi_sim = vi_double_logistic(doy_i, &curr_theta_hat.vector);

        double indv_obs = gsl_vector_get(&year_i_obs.vector, obs_i);

        sse_sim += gsl_pow_2(indv_obs - vi_sim);
      }
    }

    // 2. Draw global sigma from gamma distribution
    double sigma_tmp;
    double sigma_gamma_a, sigma_gamma_b;

    sigma_gamma_a = w->nobs / 2.0 + 0.1;
    sigma_gamma_b = sse_sim / 2.0 + 0.1;

    sigma_tmp = gsl_ran_gamma(RNG_ptr, sigma_gamma_a, 1.0 / sigma_gamma_b);

    w->sigma = 1 / sqrtl(sigma_tmp);
  }

  gsl_rng_free(RNG_ptr);
  return w;
}

/* #define gsl_vector_from_array(array, n) gsl_vector_view_array((array),
 * (n)).vector */

BLSP_Fit_T *BLSP_Fit_sample(double *obs_vec, double *doy_vec, size_t obs_count,
                      size_t *idx_vec, size_t idx_count,
                      double init_theta_mu[PARAMETER_COUNT],
                      double init_theta_sd[PARAMETER_COUNT], size_t iterations,
                      size_t burn) {

  BLSP_TimeSeries_T *X = BLSP_TimeSeries_alloc(idx_count, obs_count);
  BLSP_TimeSeries_init(X, obs_vec, doy_vec, idx_vec);

  gsl_vector_view theta_mu_view =
      gsl_vector_view_array(init_theta_mu, PARAMETER_COUNT);
  gsl_vector_view theta_sd_view =
      gsl_vector_view_array(init_theta_sd, PARAMETER_COUNT);

  BLSP_Fit_T *fit = BLSP_sampler(X, view_to_vec(theta_mu_view),
                                 view_to_vec(theta_sd_view), iterations, burn);

  BLSP_TimeSeries_free(X);
  return fit;
}

/* int BLSP_mcmc(double *obs_vec, double *doy_vec, size_t obs_count, size_t
 * *idx_vec, size_t idx_count, */
/*               double init_theta_mu[7], double init_theta_sd[7], */
/*               size_t iterations, size_t burn) { */

/*   gsl_vector v; */
/*   v = gsl_vector_view_array(obs_vec, obs_count).vector; */

/*   gsl_vector vv = gsl_vector_from_array(obs_vec, obs_count); */

/*   BLSP_TimeSeries_T *X = BLSP_TimeSeries_alloc(idx_count, obs_count); */
/*   X->time = &v; */
/*   X->time = &vv; */

/* }; */
