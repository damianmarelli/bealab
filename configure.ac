# Package config ---------------------------------------------------------------
AC_INIT([libbealab], [1.0])
AC_CONFIG_FILES([Makefile \
                 doc/Makefile \
                 include/Makefile \
                 src/Makefile \
                 ])

# Release/debug conigurations --------------------------------------------------
AC_ARG_ENABLE(debug,
    AS_HELP_STRING([--enable-debug], [Compile with debugging information.]),
    [AC_SUBST([DEBUG], [true])],
    [AC_SUBST([DEBUG], [false])])

# Checks -----------------------------------------------------------------------
AC_PREREQ([2.69])
AC_PROG_CXX

# Required Fortran compiler
AC_PROG_FC
AC_FC_SRCEXT([f],[],[echo "*** No compatible Fortran compiler found. Stop."; exit])

# Required libraries
AC_CHECK_LIB([fftw3],[main],[],[echo "*** Library fftw not found. Stop."; exit])
AC_CHECK_LIB([lapack],[main],[],[echo "*** Library lapack not found. Stop."; exit])
AC_CHECK_LIB([nlopt],[main],[],[echo "*** Library nlopt not found. Stop."; exit])
AC_CHECK_LIB([boost_serialization],[main],[],[echo "*** Library boost_serialization not found. Stop."; exit])

# Required libraries: gsl
AC_CHECK_LIB([gslcblas],[main],[],[echo "*** Library gslcblas not found. Stop."; exit])
AC_CHECK_LIB([gsl],[main],[],[echo "*** Library gsl not found. Stop."; exit])

# Required libraries: libmatio
AC_CHECK_LIB([hdf5],[main],[],[echo "*** Library hdf5 not found. Stop."; exit])
AC_CHECK_LIB([matio],[main],[],[echo "*** Library matio not found. Stop."; exit])

# Optional libraries
AC_CHECK_LIB([boost_mpi],[main],[],[echo "*** Library boost_mpi not found. MPI functionality disabled."])
AC_CHECK_LIB([ginac],[main],[],[echo "*** Library ginac not found. Symbolic module disabled."])
AC_CHECK_LIB([eng],[main],[],[echo "*** Matlab library libeng not found. Matlab interface disabled."])
AC_CHECK_LIB([GClasses],[main],[],[echo "*** Library GClasses not found. Machine learning  module disabled."])

# Optional libraries: python2.7
AC_CHECK_LIB([python2.7],[main],[],[echo "*** Library python2.7 not found. Python interface disabled."])
AC_CHECK_HEADERS([python2.7/numpy/arrayobject.h],[],[echo "*** Header python2.7/numpy/arrayobject.h not found. Python interface disabled."],[#include <python2.7/Python.h>])

# Optional libraries: libipopt
AC_CHECK_LIB([coinmetis],[main],[],[echo "*** Library coinmetis not found. Optimization via Ipopt disabled."])
AC_CHECK_LIB([coinmumps],[main],[],[echo "*** Library coinmumps not found. Optimization via Ipopt disabled."])
AC_CHECK_LIB([ipopt],[main],[],[echo "*** Library ipopt not found. Optimization via Ipopt disabled."])

# Optional external programs
AC_CHECK_PROG([GNUPLOT],[gnuplot],[-DHAVE_GNUPLOT])
AC_CHECK_PROG([MAXIMA],[maxima],[-DHAVE_MAXIMA])
AC_CHECK_PROG([PYTHON2],[python2],[-DHAVE_PYTHON2])
AC_CHECK_PROG([OCTAVE],[octave],[-DHAVE_OCTAVE])
AC_CHECK_PROG([MATLAB],[matlab],[-DHAVE_MATLAB])

# Generate the outout ----------------------------------------------------------
AC_OUTPUT
