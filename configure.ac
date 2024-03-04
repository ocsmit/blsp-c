AC_PREREQ([2.71])
AC_INIT([BLSPC], [0.0.1])
AM_INIT_AUTOMAKE([-Wall])
AM_PROG_AR
LT_INIT
AC_CONFIG_SRCDIR([src/blsp.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([m4])


dnl AC_CONFIG_MACRO_DIRS([m4])

AC_PROG_CC
AC_PROG_CXX

dnl check for required libs
dnl AC_SEARCH_LIBS([cos], [m], [], [
dnl   AC_MSG_ERROR([unable to find cos()])
dnl ])
dnl AC_SEARCH_LIBS([cblas_dgemm], [gslcblas], [], [
dnl   AC_MSG_ERROR([unable to find lib gslcblas])
dnl ])
dnl AC_SEARCH_LIBS([gsl_blas_dgemm], [gsl], [], [
dnl   AC_MSG_ERROR([unable to find lib gsl])
dnl ])

dnl AC_CHECK_HEADERS([stdint.h])

dnl AC_CHECK_HEADER([gsl/gsl_matrix.h], , [AC_MSG_ERROR(
dnl [Couldn't find the GSL header files (I searched for \
dnl <gsl/gsl_matrix.h> on the include path). If you are \
dnl using a package manager, don't forget to install the \
dnl libgsl-devel package, as well as libgsl itself.])])


AC_PATH_PROG([GSL_CONFIG], [gsl-config])
## If gsl-config was found, let's use it
if test "${GSL_CONFIG}" != ""; then
    # Use gsl-config for header and linker arguments
    GSL_CFLAGS=`${GSL_CONFIG} --cflags`
    GSL_LIBS=`${GSL_CONFIG} --libs`
else
    AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

AC_TYPE_SIZE_T

AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile src/Makefile])
#AC_CONFIG_SUBDIRS([r])

AC_CONFIG_FILES([blsp.pc])
AC_OUTPUT
