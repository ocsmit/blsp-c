AC_INIT([blspR], [0.0.1])



## Use gsl-config to find arguments for compiler and linker flags
##
## Check for non-standard programs: gsl-config(1)
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

AC_CHECK_HEADER([libblsp.h])
AC_SEARCH_LIBS([BLSP_sampler], [blsp], [BLSP_LIBS="-lblsp"], [
    AC_MSG_ERROR([BLSP_sampler not found, is blsp installed?])
])
AC_SUBST(BLSP_LIBS)


AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
