#ifndef SAMPLER_H_
#define SAMPLER_H_

#include "blsp_data.h"
gsl_matrix *blsp_sampler(BLSP_Data_t *TS_Data, double *theta_mean, double *theta_sd,
		  double sigma, size_t iterations, size_t burn);

#endif
