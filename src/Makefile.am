bin_PROGRAMS = blsp

blsp_SOURCES = test.c blsp.c timeseries.c workspace.c blsp_math.c
#LDADD = -lgsl -lgslcblas -lm
blsp_CFLAGS = $(AM_CFLAGS)

lib_LTLIBRARIES = libblsp.la
libblsp_la_SOURCES = blsp.c timeseries.c workspace.c blsp_math.c
