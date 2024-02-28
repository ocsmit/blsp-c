#ifndef BLSP_MATH_H_
#define BLSP_MATH_H_


#include "common.h"

#define M_LN_SQRT_2PI 0.918938533204672741780329736406

double vi_double_logistic(double doy, gsl_vector *parameter_vector);

double log_norm(double x, double mu, double sigma);

double log_likelihood_ratio(double candidate_value,
			    double candidate_log_likelihood, double mean_local,
			    double current_log_likelihood, double mean_global,
			    double stddev_global);


double draw_ran_gaussian(gsl_rng *RNG_ptr, double mean, double sttdev);


double vi_log_likelihood(gsl_vector *vi_vector, gsl_vector *doy_vector,
			 gsl_vector *param_vector, size_t size, double sigma);




#endif
