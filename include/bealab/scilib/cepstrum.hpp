/// @file bealab/scilib/cepstrum.hpp
/// Cepstral analysis.

#ifndef _BEALAB_CEPSTRUM_
#define	_BEALAB_CEPSTRUM_

#include <bealab/scilib/fourier.hpp>

namespace bealab
{
/// @defgroup cepstrum Cepstrum
/// Cepstral analysis.
/// @{

/// @name Exact algorithms returning a functional sequence

/// Cepstrum
template<class T>
function<complex(int)> cepstrum( const sequence<T>& x )
{
	// Functor to compute the cepstrum in the frequency domain
	function<complex(double)> ceps_f = [x]( double w ) -> complex
	{
		return log( dtft(x)( w ) );
	};

	// Return a functor that takes the idtft of the functor above
	return [ceps_f]( int t)
	{
		return idtft(ceps_f)( t );
	};
}

/// Real cepstrum
template<class T>
function<double(int)> real_cepstrum( const sequence<T>& x )
{
	// Functor to compute the real cepstrum in the frequency domain
	function<complex(double)> ceps_f = [x]( double w ) -> complex
	{
		return log( abs( dtft(x)( w ) ) );
	};

	// Return a functor that takes the idtft of the functor above
	return [ceps_f]( int t)
	{
		return idtft(ceps_f)( t );
	};
}

/// Minimum-phase cepstrum
template<class T>
function<double(int)> minimum_phase_cepstrum( const sequence<T>& x )
{
	return [x]( int t ) -> double
	{
		if( t < 0 )
			return 0;
		else if( t == 0 )
			return real_cepstrum(x)( t );
		else
			return 2 * real_cepstrum(x)( t );
	};
}

/// Inverse cepstrum
template<class T>
function<complex(int)> icepstrum( const sequence<T>& c )
{
	// Functor to compute the spectrum of the sequence
	function<complex(double)> xf = [c]( double w )
	{
		return exp( dtft(c)( w ) );
	};

	// Returm the idtft of the functor above
	return idtft(xf);
}
/// @}

/// @name DFT-based algorithms returning a sequence
/// @param N is the number of points used for DFT
/// @param a begin of the cepstrum range
/// @param b begin of the cepstrum range

/// Cepstrum
template<class T>
cseq cepstrum( const sequence<T>& x, int N, int a, int b )
{
	return idtft( log( dtft( x, N ) ) ).trunc( a, b );
}

/// Real cepstrum
template<class T>
rseq real_cepstrum( const sequence<T>& x, int N, int a, int b )
{
	return real( idtft( log( abs( dtft( x, N ) ) ) ) ).trunc( a, b );
}

/// Minimum-phase cepstrum
template<class T>
rseq minimum_phase_cepstrum( const sequence<T>& x, int N, int a, int b )
{
	rseq rceps = real_cepstrum( x, N, a, b );
	rseq mpceps = 2 * rceps;
	mpceps(0) /= 2;
	for( int t = a; t < 0; t++ )
		mpceps(t) = 0;
	return mpceps;
}

/// Inverse cepstrum
template<class T>
cseq icepstrum( const sequence<T>& c, int N, int a, int b )
{
	return idtft( exp( dtft( c, N ) ) ).trunc( a, b );
}
/// @}

/// @}
}
#endif
