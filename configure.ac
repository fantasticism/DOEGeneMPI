#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])

AC_INIT([DOEGeneMPI], [version-0.1], [joelfred@slac.stanford.edu])
AM_INIT_AUTOMAKE

PKG_CHECK_MODULES([FFTW3], [fftw3 >= 3])
AX_PROG_CXX_MPI

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
