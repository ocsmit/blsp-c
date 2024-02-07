#include "blsp_math.h"

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


double log_likelihood_ratio(double candidate_value,
			    double candidate_log_likelihood, double mean_local,
			    double current_log_likelihood, double mean_global,
			    double stddev_global) {

  // Compute values to compute accetence ratio
  double log_1 = log_norm(candidate_value, mean_global, stddev_global);
  double log_2 = log_norm(mean_local, mean_global, stddev_global);

  return (candidate_log_likelihood - current_log_likelihood) + log_1 - log_2;
}

double draw_ran_gaussian(gsl_rng *RNG_ptr, double mean, double sttdev) {
  // Draw a candidate for parameter
  return gsl_ran_gaussian(RNG_ptr, sttdev) + mean;
}

double vi_log_likelihood(gsl_vector *vi_vector, gsl_vector *doy_vector,
			 gsl_vector *param_vector, size_t size, double sigma) {
  double final_likelihood = 0, A_i, B_i;
  for (size_t i = 0; i < size; ++i) {
    A_i = gsl_vector_get(vi_vector, i);
    B_i = gsl_vector_get(doy_vector, i);
    final_likelihood += log_norm(
	gsl_vector_get(vi_vector, i),
	vi_double_logistic(gsl_vector_get(doy_vector, i), param_vector), sigma);
  }
  return final_likelihood;
}
