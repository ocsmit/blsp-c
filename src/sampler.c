#include "blsp.h"
#include "blsp_data.h"
#include "blsp_math.h"
#include "common.h"
#define PARAMETER_COUNT 7
#define BURN_STEPS 50

gsl_matrix *blsp_sampler(BLSP_Data_t *TS_Data, double *theta_mean, double *theta_sd,
		  double sigma, size_t iterations, size_t burn) {

  size_t observation_count = TS_Data->n_obs;
  size_t year_count = TS_Data->n_years;
  double m = 0.1;
  double w = year_count / (year_count + m);

  gsl_vector *theta_mean_vector;
  gsl_vector *theta_sd_vector;
  gsl_vector *current_log_like;
  gsl_vector *theta_hat_candidate_vector;

  gsl_matrix *theta_hat;
  gsl_matrix *attempts;
  gsl_matrix *acceptences;
  gsl_matrix *metropolis_hastings_matrix;


  // Create random number generator
  const gsl_rng_type *RNG_t;
  gsl_rng *RNG_ptr;
  gsl_rng_env_setup();
  RNG_t = gsl_rng_taus;
  RNG_ptr = gsl_rng_alloc(RNG_t);
  gsl_rng_set(RNG_ptr, 1);

  /* BLSP_Model_t *BLSP_Model = */
  /*     BLSP_Model_alloc(year_count, observation_count, PARAMETER_COUNT); */

  gsl_vector_view initial_theta_hat =
      gsl_vector_view_array(theta_mean, PARAMETER_COUNT);
  gsl_vector_view initial_theta_sd =
      gsl_vector_view_array(theta_sd, PARAMETER_COUNT);

  theta_mean_vector = gsl_vector_alloc(PARAMETER_COUNT);
  gsl_vector_memcpy(theta_mean_vector, &initial_theta_hat.vector);
  theta_sd_vector = gsl_vector_alloc(PARAMETER_COUNT);
  gsl_vector_memcpy(theta_sd_vector, &initial_theta_sd.vector);

  // Initialize matrix of parameter means
  // The state of each parameter for each year is stored
  theta_hat = gsl_matrix_alloc(year_count, PARAMETER_COUNT);
  for (size_t j = 0; j < year_count; ++j) {
    gsl_matrix_set_row(theta_hat, j, &initial_theta_hat.vector);
  }

  // INITIAL LOG LIKELIHOOD
  // Initialize empty (zeroed) vector to hold log likelihoods
  current_log_like = gsl_vector_calloc(year_count);
  for (size_t year = 0; year < TS_Data->n_years; ++year) {
    // Variable to store summed log likelihood for each year
    double year_i_log_likelihood = 0;
    // Get value vectors for year i

    gsl_vector_view year_i_obs = BLSP_Data_get_year_obs(TS_Data, year);
    gsl_vector_view year_i_doy = BLSP_Data_get_year_doy(TS_Data, year);

    // Get parameter vector for year i
    gsl_vector_view curr_theta_hat = gsl_matrix_row(theta_hat, year);

    double ll = vi_log_likelihood(&year_i_obs.vector, &year_i_doy.vector,
				  &curr_theta_hat.vector,
				  year_i_obs.vector.size, sigma);

    // Set log likelihood
    gsl_vector_set(current_log_like, year, year_i_log_likelihood);
  }

  // Set up matrix for attempts and for acceptences for metropolis hastings
  attempts = gsl_matrix_alloc(year_count, PARAMETER_COUNT);
  gsl_matrix_set_all(attempts, 0.1);

  acceptences = gsl_matrix_alloc(year_count, PARAMETER_COUNT);
  gsl_matrix_set_all(acceptences, 0.1);

  // Set up metropolis hastings grid and initialize
  // This stores the current standard deviation for sampling a specific
  // parameter in a specific year.
  metropolis_hastings_matrix = gsl_matrix_alloc(year_count, PARAMETER_COUNT);
  for (size_t j = 0; j < year_count; ++j) {
    gsl_matrix_set_row(metropolis_hastings_matrix, j, theta_sd_vector);
  }
  gsl_matrix_scale(metropolis_hastings_matrix, 5.0);

  // Initialize a vector to hold a /copy/ of the current year's vector so that
  // we can update the years parameter vector with the current parameters draw
  // without overwriting the value in the original theta_hat matrix until we
  // decide to accept the candidate draw.
  theta_hat_candidate_vector =
      gsl_vector_calloc(PARAMETER_COUNT); // Candidate vector

  // Declare vector views to be initalized with each year at each iteration
  gsl_vector_view year_i_obs;
  gsl_vector_view year_i_doy;
  gsl_vector_view theta_hat_year;

  for (size_t iter = 0; iter < iterations; ++iter) {
    for (size_t year = 0; year < year_count; ++year) {

      // Create vector views for the data in the current year
      year_i_obs = BLSP_Data_get_year_obs(TS_Data, year);
      year_i_doy = BLSP_Data_get_year_doy(TS_Data, year);
      theta_hat_year = gsl_matrix_row(theta_hat, year);

      for (size_t param = 0; param < PARAMETER_COUNT; ++param) {

	// Increment attempt matrix
	gsl_matrix_set(attempts, year, param,
		       gsl_matrix_get(attempts, year, param) + 1);

	// memcopy so we can change the current parameter without setting it in
	// the original
	gsl_vector_memcpy(theta_hat_candidate_vector, &theta_hat_year.vector);

	// Draw a candidate for parameter from gaussian distribution
	double candidate_draw = draw_ran_gaussian(
	    RNG_ptr, gsl_matrix_get(theta_hat, year, param),
	    gsl_matrix_get(metropolis_hastings_matrix, year, param));

	gsl_vector_set(theta_hat_candidate_vector, param, candidate_draw);

	// Compute the candidate log likelihood for current year
	double candidate_year_i_log_likelihood = vi_log_likelihood(
	    &year_i_obs.vector, &year_i_doy.vector, theta_hat_candidate_vector,
	    year_i_doy.vector.size, sigma);

	/* for (size_t obs_i = 0; obs_i < year_i_obs.vector.size; ++obs_i) { */
	/*   // Compute simulated VI */
	/*   double doy_i = gsl_vector_get(&year_i_doy.vector, obs_i); */
	/*   vi_sim = vi_double_logistic(doy_i, theta_hat_candidate_vector); */

	/*   // Draw logp(x) */
	/*   ll_draw = log_norm(gsl_vector_get(&year_i_obs.vector, obs_i),
	 * vi_sim, */
	/*		     sigma); */

	/*   // Sum up log likelihood */
	/*   candidate_year_i_log_likelihood += ll_draw; */
	/* } */

	// Compute acceptance ratio
	double R = log_likelihood_ratio(
	    candidate_draw, candidate_year_i_log_likelihood,
	    gsl_matrix_get(theta_hat, year, param),
	    gsl_vector_get(current_log_like, year),
	    gsl_vector_get(theta_mean_vector, param),
	    gsl_vector_get(theta_sd_vector, param));

	// Decide to accept or reject
	double acceptence_prob = logl(gsl_ran_flat(RNG_ptr, 0, 1));

	if (acceptence_prob < R) {
	  double current_acceptences = gsl_matrix_get(acceptences, year, param);

	  gsl_matrix_set(theta_hat, year, param, candidate_draw);
	  gsl_vector_set(current_log_like, year,
			 candidate_year_i_log_likelihood);
	  gsl_matrix_set(acceptences, year, param, current_acceptences + 1);
	}

	// Tune metropolis hastings samples
	if (iter < burn) {
	  double attempt_cnt = gsl_matrix_get(attempts, year, param);
	  double acceptence_cnt = gsl_matrix_get(acceptences, year, param);

	  // Tune every 50 iterations
	  if (attempt_cnt > BURN_STEPS) {
	    double current_sample_stddev =
		gsl_matrix_get(metropolis_hastings_matrix, year, param);

	    // If acceptence rate is low, tighten stddev
	    if (acceptence_cnt / attempt_cnt < 0.2) {
	      gsl_matrix_set(metropolis_hastings_matrix, year, param,
			     0.8 * current_sample_stddev);
	    }

	    // If acceptence rate is high, expand sttdev
	    if (acceptence_cnt / attempt_cnt > 0.6) {
	      gsl_matrix_set(metropolis_hastings_matrix, year, param,
			     1.2 * current_sample_stddev);
	    }

	    // Reset the attempt and acceptence value for current parameter
	    gsl_matrix_set(acceptences, year, param, 0);
	    gsl_matrix_set(attempts, year, param, 0);
	  }
	}
      } // End metropolis hastings for parameter vector

      // TODO: Store parameters here
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
      for (size_t year = 0; year < TS_Data->n_years; ++year) {
	double th_p = gsl_matrix_get(theta_hat, year, p);
	column_sum += th_p;
      }
      column_mean = column_sum / TS_Data->n_years;

      theta_mean_draw = draw_ran_gaussian(
	  RNG_ptr, w * column_mean + (1 - w) * theta_mean[p],
	  gsl_vector_get(theta_sd_vector, p) / sqrtl((year_count + 0.1)));

      gsl_vector_set(theta_mean_vector, p, theta_mean_draw);
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
      for (size_t year = 0; year < TS_Data->n_years; ++year) {
	xi = gsl_matrix_get(theta_hat, year, p);
	xj = gsl_vector_get(theta_mean_vector, p); // move this out (future)
	sse_i = xi - xj;
	sse_2 = gsl_pow_2(sse_i);

	stddev_sse += sse_2;
      }

      // 2. Draw new stddev from gamma dist
      double gamma_a = 0, gamma_b = 0, theta_stddev_p = 0;
      // Compute parameters
      gamma_a = year_count / 2.0 + 0.1;
      gamma_b = stddev_sse / 2.0 + 0.1;

      // Draw sample
      theta_stddev_p = gsl_ran_gamma(RNG_ptr, gamma_a, 1.0 / gamma_b);

      // intermediate variable for debugging
      double theta_stddev_p_final = 1 / sqrtl(theta_stddev_p);

      gsl_vector_set(theta_sd_vector, p, theta_stddev_p_final);
    } // END THETA STDDEV UPDATE

    // GLOBAL STDDEV **********************************************************
    // update overall standard dev (noise)
    double sse_sim = 0;
    for (size_t year = 0; year < TS_Data->n_years; ++year) {
      // 1. Calculate VI SSE for each year and add to SSE

      // Get value vectors for year i
      gsl_vector_view year_i_obs = BLSP_Data_get_year_obs(TS_Data, year);
      gsl_vector_view year_i_doy = BLSP_Data_get_year_doy(TS_Data, year);

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

    sigma_tmp = gsl_ran_gamma(RNG_ptr, sigma_gamma_a, 1 / sigma_gamma_b);

    sigma = 1 / sqrtl(sigma_tmp);
  }

  /* FILE *out_f = fopen("./data/test.dat", "wb"); */
  /* gsl_matrix_fwrite(out_f, metropolis_hastings_matrix); */
  /* fclose(out_f); */
  return theta_hat;
}
