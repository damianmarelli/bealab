/// @file bealab/extensions/signal/misc.hpp
/// Miscellaneous signal processing routines.

#ifndef _BEALAB_EXTENSIONS_SIGNAL_MISC_
#define	_BEALAB_EXTENSIONS_SIGNAL_MISC_

#include <bealab/core/blas.hpp>
#include <bealab/scilib/rootfind.hpp>
#include <bealab/extensions/control/linsys.hpp>

namespace bealab
{
namespace signal
{
//------------------------------------------------------------------------------
/// @defgroup misc Miscellaneous
/// Miscellaneous signal processing routines.
/// @{

control::state_space spectral_realization( const mat<rseq> &Rx );

/// Computes Y such that X = Y*Y.A(), within tolerance tol
mat<cseq> spectral_factorization( const mat<cseq>& X, double tol );

control::transfer_function butter( int order, double bandwidth, bool analog=false );
double raisedcosine( double t, double Fs, double beta );
double root_raisedcosine( double t, double Fs, double beta );
rseq raisedcosine( int order, double ws, double beta );
rseq root_raisedcosine( int order, double ws, double beta );
rseq hamming( int N );

/// Hertz -> Barks conversion
inline
double hertz2bark( double f )
{
	return 13 * atan(0.00076 * f) + 3.5 * pow( atan(f/7500), 2 );
};

/// Barks -> Hertz conversion
inline
double bark2hertz( double b )
{
	return fzero( [b](double f){ return hertz2bark(f) - b; }, 0, 25e3 );
}

}
}
#endif
