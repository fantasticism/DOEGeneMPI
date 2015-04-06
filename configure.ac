#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])

AC_INIT([DOEGeneMPI], [version-0.1], [joelfred@slac.stanford.edu])
AM_INIT_AUTOMAKE

PKG_CHECK_MODULES([FFTW3], [fftw3 >= 3.3.4],
		  FFTW3_LIBS="$FFTW3_LIBS -lfftw3_mpi"
		  )
PKG_CHECK_MODULES([GSL], [gsl >= 1.16])

AX_LIB_HDF5([serial])
if test "$with_hdf5" = "yes"; then
	HDF5_LIBS="$HDF5_LIBS -lhdf5_cpp"
fi

AM_PATH_PYTHON([3.4])
AX_PYTHON_MODULE([h5py])

#echo "Mymessage: $HDF5_LIBS"

AX_PROG_CXX_MPI
AC_PROG_CC

AC_CONFIG_FILES([
	Makefile
	src/Makefile
	])
AC_OUTPUT
