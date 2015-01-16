/// @file bealab/scilib/fourier.hpp
/// Discrete Fourier transform and discrete-time Fourier transform

#ifndef _BEALAB_FOURIER_
#define	_BEALAB_FOURIER_

#include <bealab/scilib/sequences/sequence.hpp>
#include <bealab/scilib/calculus.hpp>

namespace bealab
{
/// @defgroup fourier Fourier analysis
/// Discrete Fourier transform and discrete-time Fourier transform
/// @{

/// @name Time-frequency shifts
template<class T>
Seq<T> time_shift( Seq<T> X, int t )
{
    X.t1( X.t1()+t );
    return X;
}

template<class T>
Seq<decltype(T()*complex(0,0))> frequency_shift( const Seq<T>& X, double w )
{
    // Modulating function
    const rvec& t = linspace( X.t1(), X.t2(), X.size() );
    cseq A( exp(i*w*t), X.t1() );

    // Return frequency shifted version
    return element_prod( X, A );
}
/// @}

/// @name Discrete Fourier transform (DFT)

/// DFT for real vectors
Vec<complex> dft( const Vec<double>& x );

/// DFT for complex vectors
Vec<complex> dft( const Vec<complex>& x );

/// DFT for template vectors
template<class T>
Vec<decltype(T()*complex(0))> dft( const Vec<T>& X )
{
	typedef decltype(T()*complex(0)) R;
    int I = X.size();
    const rvec& W = 2*pi * vrange( 0, I ) / I;
    Vec<R> A(I);
    for( int i = 0; i < I; i++ ) {
        rvec t = vslice( 0, 1, X.size()-1 );
        cvec phase = -j*W(i)*t;
        A(i) =  sum( element_prod( X, exp(phase) ) );
    }
    return A;
}

/// Inverse DFT for real vectors
Vec<complex> idft( const Vec<double>& x );

/// Inverse DFT for complex vectors
Vec<complex> idft( const Vec<complex>& x );

/// Inverse DFT for template vectors
template<class T>
Vec<decltype(T()*complex(0))> idft( const Vec<T>& X )
{
	typedef decltype(T()*complex(0)) R;
    int I = X.size();
    const rvec& W = 2*pi * vrange( 0, I ) / I;
    Vec<R> A(I);
    for( int i = 0; i < I; i++ ) {
        const rvec& t = slice( 0, 1, X.size()-1 );
        cvec phase = j*W(i)*t;
        A(i) =  sum( element_prod( X, exp(phase) ) ) / I;
    }
    return A;
}
/// @}

/// @name Discrete-time Fourier Transform (DTFT)

/// DTFT of a sequence
template<class T>
function<decltype(T()*complex(0))(double)> dtft( const Seq<T>& x )
{
	typedef decltype(T()*complex(0)) R;
	return [x]( double w ) -> R
	{
		ivec t = vrange( x.t1(), x.t2()+1 );
		return sum( element_prod( x.vec(), exp(-j*w*t) ) );
	};
}

/// DTFT of a sequence evaluated at a vector of frequencies.
template<class T>
Vec<decltype(T()*complex(0))> dtft( const Seq<T>& X, const rvec& W )
{
	return entrywise( dtft(X), W );
}
//template<class T>
//Vec<decltype(T()*complex())> dtft( const Seq<T>& X, const rvec& W )
//{
//    int I = W.size();
//    typedef decltype(T()*complex()) R;
//    Vec<R> A(I);
//    for( int i = 0; i < I; i++ ) {
//    	auto Fs = frequency_shift( X, -W(i) );
//        A(i) = sum( Fs );
//    }
//
//    return A;
//}

/// DTFT of a sequence evaluated at a regular grid of frequencies
/// with a given number of points.
/// The grid is over [-pi,pi) and has M points.
template<class T>
Vec<decltype(T()*complex(0))> dtft( Seq<T> X, int M )
{
    int L = X.size();
    int K = ceil( (double)L/M );
    int N = K*M;

	// Deal with M = 0
    typedef decltype(T()*complex(0)) R;
    if( L == 0 )
    	return zeros(M);

    const Seq<T>& Xt   = X.trunc( X.t1(), N+X.t1()-1 );
    const Vec<T>& vec  = circshift( Xt.vec(), X.t1() );
    const Vec<R>& fvec = dft( vec );
    const ivec& idxs   = K * vrange( 0, M );
    return fvec( indirect( idxs ) );
}

/// Inverse DTFT.
template<class T>
function<decltype(T()*complex(0))(int)> idtft( const function<T(double)>& X )
{
	typedef decltype(T()*complex(0)) R;
	return [X]( int t ) -> R
	{
		function<complex(double)> integrand = [X,t](double w) {return X(w)*exp(j*w*t);};
		return 1/(2*pi) * integral( integrand, -pi, pi );
	};
}

///// Inverse DTFT specifying the interval.
//template<class T>
//cseq idtft( const function<T(double)>& X, int a, int b )
//{
//	int I = b-a+1;
//	cseq x(I,a);
//	for( int t = a; t <= b; t++ )
//		x(t) = idtft(X)(t);
//	return x;
//}
// /// Inverse DTFT.
//inline
//cseq idtft( const function<complex(double)>& X )
//{
//	// Integrand
//	bool collect;
//	rvec w_samples;
//	int t;
//	function<complex(double)> integrand = [&](double w)
//	{
//		if(collect)
//			w_samples.push_back(w);
//		return X(w)*exp(j*w*t);
//	};
//
//	// Auto-detect maximum time
//	collect = true;
//	t = 0;
//	integral( integrand, -pi, pi );
//	sort( w_samples.begin(), w_samples.end() );
//	int I = w_samples.size();
//	double delta_w = inf;
//	for( int i = 0; i < I-1; i++ ) {
//		double dw = w_samples(i+1) - w_samples(i);
//		if( dw < delta_w )
//			delta_w = dw;
//	}
//	int T = ceil(pi/delta_w);
//	(new figure)->plot(w_samples);
//	cin.get();
//
//	// Compute the idtft
//	cseq x(2*T+1,-T);
//	collect = false;
//	for( t = -T; t <= T; t++ )
//		x(t) = 1/(2*pi) * integral( integrand, -pi, pi );
//
//	return x;
//}

/// Inverse DTFT of a vector of frequency samples.
/// The vector X specifies the DFT values over a regular grid of frequencies.
/// The grid is over [-pi,pi) and has X.size() points.
/// The returned sequence is assumed to start at t1.
template<class T>
Seq<decltype(T()*complex(0))> idtft( const Vec<T>& X, int t1 )
{
    typedef decltype(T()*complex(0)) R;
    const Vec<R>& Yo = idft( X );
    const Vec<R>& Y  = circshift( Yo, -t1 );
    return Seq<R>( Y, t1 );
}

/// Inverse DTFT of a vector of frequency samples default delay.
/// The vector X specifies the DFT values over a regular grid of frequencies.
/// The grid is over [-pi,pi) and has X.size() points.
/// The returned sequence is assumed to start at -round(X.size()/2).
template<class T>
Seq<decltype(T()*complex(0))> idtft(const Vec<T>& X )
{
    int t1 = -round(X.size()/2);
    return idtft( X, t1 );
}

//------------------------------------------------------------------------------
/// @name Auxiliar functions

/// Shift zero-frequency component to the center of the spectrum
template<class T>
Vec<T> ftshift( const Vec<T>& x )
{
	int N = x.size();
	if( mod(N,2) == 0 )
		return circshift( x, N/2 );
	else
		return circshift( x, (N-1)/2 );
}

/// Inverse ftshift
template<class T>
Vec<T> iftshift( const Vec<T>& x )
{
	int N = x.size();
	if( mod(N,2) == 0 )
		return circshift( x, -N/2 );
	else
		return circshift( x, -(N-1)/2 );
}

/// Phase unwrapping
rvec unwrap( rvec p );
/// @}

//------------------------------------------------------------------------------
/// @name Discrete Fourier transform matrices

/// Direct DFT matrix
cmat dft( int N );

/// inverse DFT matrix
cmat idft( int N );
/// @}

/// @}
}
#endif
