#ifndef COMMON_H_
#define COMMON_H_


// stdlib headers
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


// GSL headers
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_vector.h>


#define PARAMETER_COUNT 7
#define BURN_STEPS 50

#define RETURN_IF_NULL(x) if (!x) { return ; }


#define view_to_vec(v) &v.vector


#endif // !COMMON_H_
