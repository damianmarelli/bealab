/// @file bealab/extensions/signal/mrsigpro.hpp
/// Routines for downsampling, upsampling and polyphase representation.

#ifndef _BEALAB_EXTENSIONS_MRSIGPRO_
#define	_BEALAB_EXTENSIONS_MRSIGPRO_

#include <bealab/scilib/sequences.hpp>

namespace bealab
{
namespace signal
{
//------------------------------------------------------------------------------
/// @defgroup mrsigpro Multirate signal processing
/// Routines for downsampling, upsampling and polyphase representation.
/// @{

/// @name Downsampling and upsampling
template <class T>
sequence<T> phase( const sequence<T> &x, int D, int d )
{
	int t1 = ceil( (x.t1()-d)/(double)D );
	int t2 = floor( (x.t2()-d)/(double)D );
	int I  = t2-t1+1;
	vec<T> v(I);
	for( int i = 0; i < I; i++ )
		v(i) = x( (i+t1)*D + d );
	return sequence<T>( v, t1 );
}

template <class T>
sequence<T> downsample( const sequence<T> &x, int D )
{
	return phase( x, D, 0 );
}

template <class T>
vec<sequence<T>> downsample( const vec<sequence<T>> &x, int D )
{
	return x.apply( [D](const rseq &s){return downsample(s,D);} );
}

template <class T>
mat<sequence<T>> downsample( const mat<sequence<T>> &x, int D )
{
	return x.apply( [D](const rseq &s){return downsample(s,D);} );
}

template <class T>
sequence<T> upsample( const sequence<T> &x, int D )
{
	if( x.size() == 0 )
		return sequence<T>();
	int xt1    = x.t1();
	int xt2    = x.t2();
	int xI     = xt2-xt1+1;
	int t1     = D*xt1;
	int t2     = D*xt2;
	int I      = t2-t1+1;
	vec<T> v   = zeros(I);
	const vec<T> &x_ = x.buffer();
	for( int i = 0; i < xI; i++ )
		v(i*D) = x_(i);
	return sequence<T>( v, t1 );
}

template <class T>
vec<sequence<T>> upsample( const vec<sequence<T>> &x, int D )
{
	return x.apply( [D](const rseq &s){return upsample(s,D);} );
}

template <class T>
mat<sequence<T>> upsample( const mat<sequence<T>> &x, int D )
{
	return x.apply( [D](const rseq &s){return upsample(s,D);} );
}
/// @}

/// @name Polyphase representation
template <class T>
vec<sequence<T>> fb2pp_signal( const sequence<T> &x, int D )
{
	vec<sequence<T>> X(D);
	for( int d = 0; d < D; d++ )
		X(d) = phase( x, D, -d );
	return X;
}

template <class T>
sequence<T> pp2fb_signal( const vec<sequence<T>> &X )
{
	int D = X.size();
	sequence<T> x;
	for( int d = 0; d < D; d++ )
		x += time_shift( upsample( X(d), D ), -d );
	return x;
}

template <class T>
mat<sequence<T>> fb2pp_filterbank( const vec<sequence<T>> &h, int D )
{
	int M = h.size();
	mat<sequence<T>> H(M,D);
	for( int m = 0; m < M; m++ )
		for( int d = 0; d < D; d++ )
			H(m,d) = phase( h(m), D, d );
	return H;
}

template <class T>
vec<sequence<T>> pp2fb_filterbank( const mat<sequence<T>> &H )
{
	int M = H.rows();
	int D = H.cols();
	vec<sequence<T>> h(M);
	for( int m = 0; m < M; m++ )
		for( int d = 0; d < D; d++ )
			h(m) += time_shift( upsample( H(m,d), D ), d );
	return h;
}

template <class T>
mat<sequence<T>> fb2pp_system( const sequence<T> &g, int D )
{
	mat<sequence<T>> G(D,D);
	for( int d = 0; d < D; d++ )
		for( int e = 0; e < D; e++ )
			G(d,e) = phase( g, D, e-d );
	return G;
}

template <class T>
mat<sequence<T>> fb2pp_system( const vec<sequence<T>> &g )
{
	int D = g.size();
	mat<sequence<T>> G(D,D);
	for( int d = 0; d < D; d++ )
		for( int e = 0; e < D; e++ )
			G(d,e) = phase( g(d), D, e-d );
	return G;
}

template <class T>
vec<sequence<T>> pp2fb_system( const mat<sequence<T>> &G )
{
	if( G.size1() != G.size2() )
		error("pp2fb_system() - The polyphase matrix needs to be square");
	int D = G.size1();
	vec<sequence<T>> g(D);
	for( int d = 0; d < D; d++ ) {
		for( int e = 0; e < D; e++ ) {
			const sequence<T>& ph = upsample( G(d,e), D );
			g(d) += time_shift( ph, e-d );
		}
	}
	return g;
}
/// @}

/// @}
}
}
#endif
