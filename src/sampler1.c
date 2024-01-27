//#define DEBUG
#include "common.h"
#include "ts_data.h"
#include <stdio.h>

#define GSL_RNG_SEED 1
#define W 0.9950249 //     w <- nyear / (nyear + m) # coeff for post norm

#define M_LN_SQRT_2PI 0.918938533204672741780329736406
#define N_PARAMETER 7

// GENERAL FUNKS **************************************************************

/*
 * Double logistic function with greendown parameter (Elmore et al., 2012)
 */
double vi_double_logistic(double doy, gsl_vector *parameter_vector) {

  double m1_l = gsl_cdf_logistic_P(gsl_vector_get(parameter_vector, 0), 1);
  double m7_l = gsl_cdf_logistic_P(gsl_vector_get(parameter_vector, 6), 1);

  double term1, // (m2 - m7 * t)
      term2,    // 1 / (1 + e^((m3 - t) / m4))
      term3;    // 1 / (1 + e^((m5 - t) / m6))
  int status = 0;

  gsl_sf_result_e10 term2_exp;
  gsl_sf_result_e10 term3_exp;
  term1 = gsl_vector_get(parameter_vector, 1) - m7_l * doy;

  gsl_sf_exp_e10_e((gsl_vector_get(parameter_vector, 2) - doy) /
                       gsl_vector_get(parameter_vector, 3),
                   &term2_exp);
  term2 = (1 / (1 + term2_exp.val));

  gsl_sf_exp_e10_e((gsl_vector_get(parameter_vector, 4) - doy) /
                       gsl_vector_get(parameter_vector, 5),
                   &term3_exp);

  term3 = (1 / (1 + term3_exp.val));

  return m1_l + term1 * (term2 - term3);
}

/*
 * Safely compute /log/ probabilty from normal distribution.
 * Adapated from dnorm.c & dnorm4() from R Core.
 */
double log_norm(double x, double mu, double sigma) {
  double xx;
  xx = (x - mu) / sigma;
  xx = fabs(xx);
  return -(M_LN_SQRT_2PI + 0.5 * xx * xx + gsl_sf_log(sigma));
}

// MAIN ***********************************************************************
int main(int argc, char const *argv[]) {

  // Parameters taken from simulated data in R
  double theta_mean[N_PARAMETER] = {-1.998392,  0.960355, 120.702350, 9.263498,
                                    288.853856, 9.166685, -6.592421};
  double theta_sd[N_PARAMETER] = {0.07057906, 0.05609551, 1.08944966,
                                  0.88183154, 1.55979462, 1.20727157,
                                  0.19881890};

  // Copy over theta mean and theta standard deviation to GSL_Vectors
  // Could also just create a GSL array view to the two arrays but would rather
  // just copy the data for now
  gsl_vector *theta_sd_view = gsl_vector_calloc(N_PARAMETER);
  gsl_vector *theta_mean_view = gsl_vector_calloc(N_PARAMETER);
  for (size_t i = 0; i < N_PARAMETER; ++i) {
    gsl_vector_set(theta_mean_view, i, theta_mean[i]);
    gsl_vector_set(theta_sd_view, i, theta_sd[i]);
  }

  // Set up RNG
  const gsl_rng_type *RNG_t;
  gsl_rng *RNG_ptr;
  // Create random number generator
  gsl_rng_env_setup();
  RNG_t = gsl_rng_taus;
  RNG_ptr = gsl_rng_alloc(RNG_t);
  gsl_rng_set(RNG_ptr, 1);

  // We simulated 20 years and 334 obs
  size_t num_years = 20, num_obs = 334, curr_yr_id = 0;
  double prior_strength = 0.1; // used when updating MCMC

  // GSL vectors to store data
  gsl_vector *doy_vector = gsl_vector_alloc(num_obs);
  gsl_vector *y_vector = gsl_vector_alloc(num_obs);
  gsl_vector *yr_id_vector = gsl_vector_alloc(num_years);

  plog("Reading data");
  // Read in test data. GLS will throw error if file not found
  FILE *file_doy = fopen("./data/doy_full.dat", "r");
  gsl_vector_fscanf(file_doy, doy_vector);
  fclose(file_doy);

  FILE *file_y = fopen("./data/y_full.dat", "r");
  gsl_vector_fscanf(file_y, y_vector);
  fclose(file_y);

  FILE *file_yr_id = fopen("./data/year_idx.dat", "r");
  gsl_vector_fscanf(file_y, yr_id_vector);
  fclose(file_yr_id);

  // Store our initial state. This is used in the sampler but *does not change*
  gsl_vector *initial_theta_hat = gsl_vector_alloc(N_PARAMETER);
  gsl_vector_memcpy(initial_theta_hat, theta_mean_view);

  // Allocate structure to hold data allowing us to slice data by year
  TS_Data_t *TS_Data = TS_Data_alloc(num_years, num_obs);

  // Change pointer to the needed vectors
  TS_Data->obs = y_vector;
  TS_Data->doy = doy_vector;

  // Set up the year_idx vector. Important to note that we iterate from one so
  // that we make sure the year indexing values start from 0.
  // See ts_data.h for more info
  plog("Set up year idx");
  for (size_t i = 1; i < TS_Data->n_years +1; ++i) {
    gsl_vector_set(TS_Data->year_idx, i,
                   gsl_vector_get(TS_Data->year_idx, i - 1) +
                       gsl_vector_get(yr_id_vector, i - 1));
  }



  // Initialize matrix of parameter means.
  // The parameter vector state for each year is updated in this matrix,
  // and the global mean and stddev are updated using the values stored here
  // as well.
  plog("initializing theta hat");
  gsl_matrix *theta_hat = gsl_matrix_alloc(num_years, N_PARAMETER);
  for (size_t j = 0; j < num_years; ++j) {
    gsl_matrix_set_row(theta_hat, j, theta_mean_view);
  }

  double sigma = 0.08040725;

  // INITIAL LOG LIKELIHOOD ***************************************************
  // Initial Log likelihood at the current value of theta for each year.

  // Initialize empty (zeroed) vector to hold log likelihoods
  gsl_vector *current_log_like = gsl_vector_calloc(num_years);

  // Iterate over each year

  plog("Computing initial log likelihoods");
  for (size_t year = 0; year < num_years; ++year) {
    // Variable to store summed log likelihood for each year
    double year_i_log_likelihood = 0;
    // Get value vectors for year i

    gsl_vector_view year_i_obs = TS_Data_get_year_obs(TS_Data, year);
    gsl_vector_view year_i_doy = TS_Data_get_year_doy(TS_Data, year);

    // Get parameter vector for year i
    gsl_vector_view curr_theta_hat = gsl_matrix_row(theta_hat, year);

    // Iterate over each obs in year
    double vi_sim = 0;
    double ll_draw = 0;

    for (size_t obs_i = 0; obs_i < year_i_obs.vector.size; ++obs_i) {
      // Compute simulated VI
      double doy_i = gsl_vector_get(&year_i_doy.vector, obs_i);
      vi_sim = vi_double_logistic(doy_i, &curr_theta_hat.vector);

      // Draw logp(x)
      double vi_true = gsl_vector_get(&year_i_obs.vector, obs_i);
      ll_draw = log_norm(vi_true, vi_sim, sigma);

      // Sum up log likelihood
      year_i_log_likelihood += ll_draw;
    }

    // Set log likelihood
    gsl_vector_set(current_log_like, year, year_i_log_likelihood);
  }

  plog("Initializing metropolis hastings acceptence and attempt grids");
  // Set up matrix for attempts and for acceptences for metropolis hastings
  gsl_matrix *attempts = gsl_matrix_alloc(num_years, N_PARAMETER);
  gsl_matrix_set_all(attempts, 0.1);

  gsl_matrix *acceptences = gsl_matrix_alloc(num_years, N_PARAMETER);
  gsl_matrix_set_all(acceptences, 0.1);

  // Set up metropolis hastings grid and initialize
  // This stores the current standard deviation for sampling a specific
  // parameter in a specific year.
  plog("Initializing metropolis hastings sampling grid");
  gsl_matrix *mh_matrix = gsl_matrix_alloc(num_years, N_PARAMETER);
  for (size_t j = 0; j < num_years; ++j) {
    gsl_matrix_set_row(mh_matrix, j, theta_sd_view);
  }
  gsl_matrix_scale(mh_matrix, 5.0);

  unsigned int iterations = 7000;
  unsigned int burn = 2000;

  // NOTE: Right now we are not keeping track each parameter draw.
  //       I am diagnosing the sampler based off the final SD and mean grids.
  // gsl_matrix *parameter_tracker =
  //     gsl_matrix_calloc(num_years * 7000, N_PARAMETER);
  // gsl_matrix *parameter_tracker_1 = gsl_matrix_calloc(num_years,
  // N_PARAMETER);

  // Initialize a vector to hold a /copy/ of the current year's vector so that
  // we can update the years parameter vector with the current parameters draw
  // without overwriting the value in the original theta_hat matrix until we
  // decide to accept the candidate draw.
  gsl_vector *candidate_theta_hat_vector =
      gsl_vector_calloc(N_PARAMETER); // Candidate vector

  // **** BEGIN ***************************************************************
  plog("Starting sampler");
  for (size_t iter = 0; iter < iterations; ++iter) {
    for (size_t year = 0; year < TS_Data->n_years; ++year) {
      // Create vector views for the data in the current year
      gsl_vector_view year_i_obs = TS_Data_get_year_obs(TS_Data, year);
      gsl_vector_view year_i_doy = TS_Data_get_year_doy(TS_Data, year);
      gsl_vector_view theta_hat_year = gsl_matrix_row(theta_hat, year);

      for (size_t param = 0; param < N_PARAMETER; ++param) {

        plog("ITER: %lu\t | YEAR: %lu\t | PARAM: %lu", iter, year, param);

        // attempts ++
        gsl_matrix_set(attempts, year, param,
                       gsl_matrix_get(attempts, year, param) + 1);

        // memcopy so we can change the current parameter without setting it in
        // the original
        // Get current years parameter vector
        // gsl_vector_view candidate_theta = gsl_matrix_row(theta_hat, year);
        // gsl_vector_memcpy(candidate_vector, &candidate_theta.vector);
        // gsl_matrix_get_row(candidate_theta_hat_vector, theta_hat, year);
        gsl_vector_memcpy(candidate_theta_hat_vector, &theta_hat_year.vector);

        double cc_sd = gsl_matrix_get(mh_matrix, year, param); // Current stddev
        double cc_mn = gsl_matrix_get(theta_hat, year, param); // current mean

        // Draw a candidate for parameter
        double candidate_draw = gsl_ran_gaussian(RNG_ptr, cc_sd) + cc_mn;
        gsl_vector_set(candidate_theta_hat_vector, param, candidate_draw);

        // Compute the candidate log likelihood for year
        double vi_sim = 0;
        double ll_draw = 0;
        double candidate_year_i_log_likelihood = 0;
        for (size_t obs_i = 0; obs_i < year_i_obs.vector.size; ++obs_i) {
          // Compute simulated VI
          double doy_i = gsl_vector_get(&year_i_doy.vector, obs_i);
          vi_sim = vi_double_logistic(doy_i, candidate_theta_hat_vector);

          // Draw logp(x)
          double vi_true = gsl_vector_get(&year_i_obs.vector, obs_i);
          ll_draw = log_norm(vi_true, vi_sim, sigma);

          // Sum up log likelihood
          candidate_year_i_log_likelihood += ll_draw;
        }

        // Get the global mean and stddev for current parameter
        double theta_mn_j, theta_sd_j;
        theta_mn_j = gsl_vector_get(theta_mean_view, param);
        theta_sd_j = gsl_vector_get(theta_sd_view, param);

        // Compute values to compute accetence ratio
        double theta_i = gsl_matrix_get(theta_hat, year, param);
        double log_1 = log_norm(candidate_draw, theta_mn_j, theta_sd_j);
        double log_2 = log_norm(theta_i, theta_mn_j, theta_sd_j);

        double cur_ll = gsl_vector_get(current_log_like, year);
        double ll_dif = candidate_year_i_log_likelihood - cur_ll;

        // Compute acceptance ratio
        double R = ll_dif + log_1 - log_2;

        // Decide to accept or reject
        double ttt = gsl_ran_flat(RNG_ptr, 0, 1); // TODO: change shit variable
                                                  //       name
        double l_ttt = logl(ttt);
        // double ttt = logf((rand()/(double)(RAND_MAX)) * 2 - 1);

        // int status = gsl_fcmp(R, gsl_sf_log(gsl_ran_flat(r, 0, 1)),
        // 0.000001f);
        // if (status > 0) {
        if (l_ttt < R) { // Update values for current parameter
          double c_a = gsl_matrix_get(acceptences, year, param);

          gsl_matrix_set(theta_hat, year, param, candidate_draw);
          gsl_vector_set(current_log_like, year,
                         candidate_year_i_log_likelihood);
          gsl_matrix_set(acceptences, year, param, c_a + 1);
        }

        // Tune metropolis hastings sampler
        if (iter < burn) {
          double attempt_cnt = gsl_matrix_get(attempts, year, param);
          double acceptence_cnt = gsl_matrix_get(acceptences, year, param);

          // Tune every 50 iterations
          if (attempt_cnt > 50) {
            double curr_mh = gsl_matrix_get(mh_matrix, year, param);
            // NOTE: mh_matrix is not converging correctly

            // If acceptence rate is low, tighten stddev
            if (acceptence_cnt / attempt_cnt < 0.2) {
              gsl_matrix_set(mh_matrix, year, param, 0.8 * curr_mh);
            }

            // If acceptence rate is high, expand sttdev
            if (acceptence_cnt / attempt_cnt > 0.6) {
              gsl_matrix_set(mh_matrix, year, param, 1.2 * curr_mh);
            }

            // Reset the attempt and acceptence value for current parameter
            gsl_matrix_set(acceptences, year, param, 0);
            gsl_matrix_set(attempts, year, param, 0);
          }
        }
      } // END metropolis hastings for parameter vector

      // NOTE: When storing parameters, here is where to do it

    }   // END year loop

    // **** GIBBS *************************************************************
    // Sample for global mean and stddev

    // THETA MEAN *************************************************************
    // Update the vector storing global mean for each parameter
    // theta_mn <- rnorm(npar, w * colMeans(theta_hat) + (1 - w) * mu,
    // theta_sd / sqrt((nyear + m)))
    //
    // NOTE: When refactoring, see if it's faster to use BLAS to calculate
    // column means instead of iterating fully over each column:
    // gsl_blas_dgemv(CblasTrans, 1.0 / theta_hat->size1, theta_hat,
    // scale_vector, 0.0, col_means);

    for (size_t p = 0; p < N_PARAMETER; ++p) {
      double column_mean = 0, column_sum = 0;
      double theta_mn_mu = 0, theta_mn_sig = 0;
      double theta_mean_p = 0;

      // calculate column mean for parameter p
      for (size_t year = 0; year < TS_Data->n_years; ++year) {
        column_sum += gsl_matrix_get(theta_hat, year, p);
      }
      column_mean = column_sum / TS_Data->n_years;

      // calculate mu parameter for distribution
      theta_mn_mu =
          W * column_mean + (1 - W) * gsl_vector_get(initial_theta_hat, p);

      // calculate sigma parameter for distribution
      theta_mn_sig =
          gsl_vector_get(theta_sd_view, p) / sqrtl((num_years / 0.1));

      // Pull draw
      theta_mean_p = gsl_ran_gaussian(RNG_ptr, theta_mn_sig);
      // Apply mean transformation
      theta_mean_p += theta_mn_mu;

      // Set theta_mn[p]
      gsl_vector_set(theta_mean_view, p, theta_mean_p);
    } // END THETA MEAN UPDATE

    // THETA STDDEV ***********************************************************
    // Update the vector storing global stddev for each parameter
    // 1. Calculate SSE
    //  sse_theta_sd <- colSums(Rfast::eachrow(theta_hat, theta_mn, "-")^2)
    //
    // 2. Draw new stddev from gamma dist
    // theta_sd <- 1 / sqrt(rgamma(npar, nyear / 2 + a, sse_theta_sd / 2 + b))
    // message(theta_sd)

    for (size_t p = 0; p < N_PARAMETER; ++p) {
      // 1. SSE for parameter p
      double xi = 0, xj = 0, stddev_sse = 0;
      for (size_t year = 0; year < TS_Data->n_years; ++year) {
        xi = gsl_matrix_get(theta_hat, year, p);
        xj = gsl_vector_get(theta_mean_view, p); // move this out (future)
        double sse_i = xi - xj;
        stddev_sse += sse_i * sse_i;
      }

      // 2. Draw new stddev from gamma dist
      double gamma_a = 0, gamma_b = 0, theta_stddev_p = 0;
      // Compute parameters
      gamma_a = num_years / 2.0 + 0.1;
      gamma_b = stddev_sse / 2.0 + 0.1;

      // Draw sample
      theta_stddev_p = gsl_ran_gamma(RNG_ptr, gamma_a, gamma_b);

      // intermediate variable for debugging
      double theta_stddev_p_final = 1 / sqrtl(theta_stddev_p);

      gsl_vector_set(theta_sd_view, p, theta_stddev_p_final);
    } // END THETA STDDEV UPDATE

    // GLOBAL STDDEV **********************************************************
    // update overall standard dev (noise)
    //
    // this line just splits the theta_hat matrix into parameter vectors
    //    t1 <- lapply(seq_len(nyear), function(i) theta_hat[i, ])
    // this line computes the VI updates
    //    t2 <- unlist(mapply(f, t, t1), F, F)
    // this line computes the total SSE
    //    sse_sig <- sum((y_long - t2)^2)
    //    sigma <- 1 / sqrt(rgamma(1, (tobs) / 2 + a, sse_sig / 2 + b))
    double sse_sim = 0;
    for (size_t year = 0; year < TS_Data->n_years; ++year) {
      // 1. Calculate VI SSE for each year and add to SSE

      // Get value vectors for year i
      gsl_vector_view year_i_obs = TS_Data_get_year_obs(TS_Data, year);
      gsl_vector_view year_i_doy = TS_Data_get_year_doy(TS_Data, year);

      // We use the scale vector as a temporary place to memcpy our parameter
      // vector for year i
      // gsl_matrix_get_row(candidate_vector, theta_hat, year);

      gsl_vector_view curr_theta_hat = gsl_matrix_row(theta_hat, year);
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

    sigma_gamma_a = TS_Data->n_obs / 2.0 + 0.1;
    sigma_gamma_b = sse_sim / 2.0 + 0.1;

    sigma_tmp = gsl_ran_gamma(RNG_ptr, sigma_gamma_a, sigma_gamma_b);

    sigma = 1 / sqrtl(sigma_tmp);

  } // END iter loop

  FILE *out_f = fopen("./data/test.dat", "wb");
  // gsl_matrix_fprintf(out_f, theta_hat, "%f");
  gsl_matrix_fwrite(out_f, mh_matrix);
  fclose(out_f);


  TS_Data_free(TS_Data);
  // these are freed by TS_Data_free currently
  // gsl_vector_free(doy_vector);
  // gsl_vector_free(y_vector);
  // gsl_vector_free(yr_id_vector);

  gsl_vector_free(theta_mean_view);
  gsl_vector_free(theta_sd_view);
  gsl_vector_free(initial_theta_hat);
  gsl_vector_free(current_log_like);
  gsl_vector_free(candidate_theta_hat_vector);

  gsl_matrix_free(theta_hat);
  gsl_matrix_free(attempts);
  gsl_matrix_free(acceptences);
  gsl_matrix_free(mh_matrix);


  return 0;
}
