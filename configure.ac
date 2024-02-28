AC_PREREQ([2.71])
AC_INIT([BLSPC], [0.0.1])
AM_INIT_AUTOMAKE([-Wall -Werror])
AC_CONFIG_SRCDIR([src/blsp.c])
AC_CONFIG_HEADERS([config.h])

AC_PROG_CC
AC_PROG_CXX

dnl check for required libs
AC_SEARCH_LIBS([cos], [m], [], [
  AC_MSG_ERROR([unable to find cos()])
])
AC_SEARCH_LIBS([cblas_dgemm], [gslcblas], [], [
  AC_MSG_ERROR([unable to find lib gslcblas])
])
AC_SEARCH_LIBS([gsl_blas_dgemm], [gsl], [], [
  AC_MSG_ERROR([unable to find lib gsl])
])

AC_CHECK_HEADERS([stdint.h])

AC_CHECK_HEADER([gsl/gsl_matrix.h], , [AC_MSG_ERROR(
[Couldn't find the GSL header files (I searched for \
<gsl/gsl_matrix.h> on the include path). If you are \
using a package manager, don't forget to install the \
libgsl-devel package, as well as libgsl itself.])])

AC_TYPE_SIZE_T

AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile src/Makefile])
AC_CONFIG_SUBDIRS([r])

AC_OUTPUT
