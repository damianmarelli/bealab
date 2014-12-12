/// @file bealab/scilib/rootfind.hpp
/// Root-finding routines.

#ifndef _BEALAB_ROOTFIND_
#define	_BEALAB_ROOTFIND_

#include <bealab/core/blas.hpp>

namespace bealab
{
/// @defgroup rootfind Root-finding
/// Root-finding routines.
/// @{

/// Finds the zero of a function within a given interval.
double fzero( function<double(double)> fun, double x_lo, double x_hi,
		int max_iter=1000, double epsabs = 1e-6, double epsrel = 1e-6 );

/// @}
}
#endif
