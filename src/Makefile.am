check_PROGRAMS  = current
# gsltest_SOURCES = gsl-testing.cpp
# gsltest_LDADD   = $(GSL_LIBS)

current_SOURCES = work.cpp
current_LDADD   = $(GSL_LIBS) $(FFTW3_LIBS) $(HDF5_LIBS) $(HDF5_LDFLAGS)
current_CPPFLAGS = $(HDF5_CPPFLAGS)

clean-local:
	-rm -rf *.h5

# TESTS           = $(check_PROGRAMS)
