GSL_CFLAGS = @GSL_CFLAGS@
GSL_LIBS   = @GSL_LIBS@

bin_PROGRAMS = blsp

blsp_SOURCES = test.c blsp.c timeseries.c workspace.c blsp_math.c
blsp_CFLAGS  = $(AM_CFLAGS) $(GSL_CFLAGS)
blsp_LDFLAGS = $(AM_LDFLAGS) $(GSL_LIBS)



lib_LTLIBRARIES    = libblsp.la

libblsp_la_SOURCES = blsp.c timeseries.c workspace.c blsp_math.c
libblsp_la_CFLAGS  = $(AM_CFLAGS) $(GSL_CFLAGS)
libblsp_la_LDFLAGS = $(AM_LDFLAGS) $(GSL_LIBS)
