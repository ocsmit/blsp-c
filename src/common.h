#ifndef COMMON_H_
#define COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

// stdlib headers
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
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

// General logger
#ifdef DEBUG
#define plog(format, ...) ( {  \
  fprintf(stderr, "\033[1;34m[plog]\033[0m " format " \033[2m(%s %d)\033[0m\n", \
    ##__VA_ARGS__, __FUNCTION__, __LINE__); \
  fflush(stderr); \
})
#define perrlog(format, ...) (({ \
  fprintf(stderr, "error: " format " (%s:%d)\n", \
    ##__VA_ARGS__, __FILE__, __LINE__); \
  fflush(stderr); \
}))
#else
#define plog(format, ...)
#define perrlog(format, ...)
#endif

#ifdef __cplusplus
}
#endif


#endif // !COMMON_H_
