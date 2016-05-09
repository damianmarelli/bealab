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
AC_LANG_PUSH([C++])

# Required Fortran compiler
AC_PROG_FC
AC_FC_SRCEXT([f],[],[AC_MSG_ERROR(No compatible Fortran compiler found.)])

# Required libraries
AC_CHECK_LIB([fftw3],[main],[],[AC_MSG_ERROR(Library fftw not found.)])
AC_CHECK_LIB([lapack],[main],[],[AC_MSG_ERROR(Library lapack not found.)])
AC_CHECK_LIB([nlopt],[main],[],[AC_MSG_ERROR(Library nlopt not found.)])
AC_CHECK_LIB([gsl],[main],[],[AC_MSG_ERROR(Library gsl not found.)])
AC_CHECK_LIB([matio],[main],[],[AC_MSG_ERROR(Library matio not found.)])

# Required libraries: boost
AC_CHECK_LIB([boost_serialization],[main],[],[AC_MSG_ERROR(Library boost_serialization not found.)])
AC_CHECK_HEADER([boost/math/distributions/normal.hpp],[],[AC_MSG_ERROR(Header boost/math/distributions/normal.hpp not found.)])

# Optional libraries
#AC_CHECK_LIB([boost_mpi],[main],[],[AC_MSG_WARN(Library boost_mpi not found. MPI functionality disabled.)])
AC_CHECK_LIB([ginac],[main],[],[AC_MSG_WARN(Library ginac not found. Symbolic module disabled.)])
AC_CHECK_LIB([GClasses],[main],[],[AC_MSG_WARN(Library GClasses not found. Machine learning  module disabled.)])
AC_CHECK_LIB([ipopt],[main],[],[AC_MSG_WARN(Library ipopt not found. Optimization via Ipopt disabled.)])

# Optional libraries: python2.7
AC_CHECK_LIB([python2.7],[main],[],[AC_MSG_WARN(Library python2.7 not found. Python interface disabled.)])
AC_CHECK_HEADERS([python2.7/Python.h],[],[AC_MSG_WARN(Header python2.7/Python.h not found. Python interface disabled.)])

# Optional external programs
AC_CHECK_PROG([GNUPLOT],[gnuplot],[-DHAVE_GNUPLOT])
AC_CHECK_PROG([MAXIMA],[maxima],[-DHAVE_MAXIMA])
AC_CHECK_PROG([PYTHON2],[python2],[-DHAVE_PYTHON2])
AC_CHECK_PROG([OCTAVE],[octave],[-DHAVE_OCTAVE])
AC_CHECK_PROG([MATLAB],[matlab],[-DHAVE_MATLAB])

# Generate the outout ----------------------------------------------------------
AC_OUTPUT
