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

template<class VAL>
Mat<VAL> convmat(const Vec<VAL> &X, int J)
{
    int L = X.size();
    int I = L + J - 1;
    Mat<VAL> M = zeros(I,J);
    Vec<VAL> Z = zeros(J-1);
    Vec<VAL> B = { Z, flip(X), Z };

    for(int i=0; i<I; i++)
        M.row(i) = B.range( L+J-2-i, L+2*J-3-i );

    return M;
}

control::state_space spectral_realization( const Mat<rseq> &Rx );

/// Computes Y such that X = Y*Y.A(), within tolerance tol
Mat<cseq> spectral_factorization( const Mat<cseq>& X, double tol );

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

///// Compute the N principal components of the vectors in the columns of the
///// matrix X.
//template<class T>
//Mat<T> pca( const Mat<T>& X, int N )
//{
//	// Compute principal components
//	Mat<T> XX = X * trans(X);
//	auto DU   = eig( XX );
//	rvec d    = real(diag( get<0>(DU) ));
//	rmat U    = real(get<1>(DU));
//	int I     = X.size1();
//	Mat<T> pcs(I,N);
//	for( int n = 0; n < N; n++ )
//		pcs.column(n) = U.column(n);
//	return pcs;
//}
/// @}
}
}
#endif
