// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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
Seq<T> phase( const Seq<T> &x, int D, int d )
{
	int t1 = ceil( (x.t1()-d)/(double)D );
	int t2 = floor( (x.t2()-d)/(double)D );
	int I  = t2-t1+1;
	Vec<T> v(I);
	for( int i = 0; i < I; i++ )
		v(i) = x( (i+t1)*D + d );
	return Seq<T>( v, t1 );
}

template <class T>
Seq<T> downsample( const Seq<T> &x, int D )
{
	return phase( x, D, 0 );
}

template <class T>
Vec<Seq<T>> downsample( const Vec<Seq<T>> &x, int D )
{
	return x.apply( [D](const rseq &s){return downsample(s,D);} );
}

template <class T>
Mat<Seq<T>> downsample( const Mat<Seq<T>> &x, int D )
{
	return x.apply( [D](const rseq &s){return downsample(s,D);} );
}

template <class T>
Seq<T> upsample( const Seq<T> &x, int D )
{
	if( x.size() == 0 )
		return Seq<T>();
	int xt1    = x.t1();
	int xt2    = x.t2();
	int xI     = xt2-xt1+1;
	int t1     = D*xt1;
	int t2     = D*xt2;
	int I      = t2-t1+1;
	Vec<T> v   = zeros(I);
	const Vec<T> &x_ = x.vec();
	for( int i = 0; i < xI; i++ )
		v(i*D) = x_(i);
	return Seq<T>( v, t1 );
}

template <class T>
Vec<Seq<T>> upsample( const Vec<Seq<T>> &x, int D )
{
	return x.apply( [D](const rseq &s){return upsample(s,D);} );
}

template <class T>
Mat<Seq<T>> upsample( const Mat<Seq<T>> &x, int D )
{
	return x.apply( [D](const rseq &s){return upsample(s,D);} );
}
/// @}

/// @name Polyphase representation
template <class T>
Vec<Seq<T>> fb2pp_signal( const Seq<T> &x, int D )
{
	Vec<Seq<T>> X(D);
	for( int d = 0; d < D; d++ )
		X(d) = phase( x, D, -d );
	return X;
}

template <class T>
Seq<T> pp2fb_signal( const Vec<Seq<T>> &X )
{
	int D = X.size();
	Seq<T> x;
	for( int d = 0; d < D; d++ )
		x += time_shift( upsample( X(d), D ), -d );
	return x;
}

template <class T>
Mat<Seq<T>> fb2pp_filterbank( const Vec<Seq<T>> &h, int D )
{
	int M = h.size();
	Mat<Seq<T>> H(M,D);
	for( int m = 0; m < M; m++ )
		for( int d = 0; d < D; d++ )
			H(m,d) = phase( h(m), D, d );
	return H;
}

template <class T>
Vec<Seq<T>> pp2fb_filterbank( const Mat<Seq<T>> &H )
{
	int M = H.rows();
	int D = H.cols();
	Vec<Seq<T>> h(M);
	for( int m = 0; m < M; m++ )
		for( int d = 0; d < D; d++ )
			h(m) += time_shift( upsample( H(m,d), D ), d );
	return h;
}

template <class T>
Mat<Seq<T>> fb2pp_system( const Seq<T> &g, int D )
{
	Mat<Seq<T>> G(D,D);
	for( int d = 0; d < D; d++ )
		for( int e = 0; e < D; e++ )
			G(d,e) = phase( g, D, e-d );
	return G;
}

template <class T>
Mat<Seq<T>> fb2pp_system( const Vec<Seq<T>> &g )
{
	int D = g.size();
	Mat<Seq<T>> G(D,D);
	for( int d = 0; d < D; d++ )
		for( int e = 0; e < D; e++ )
			G(d,e) = phase( g(d), D, e-d );
	return G;
}

template <class T>
Vec<Seq<T>> pp2fb_system( const Mat<Seq<T>> &G )
{
	if( G.size1() != G.size2() )
		error("pp2fb_system() - The polyphase matrix needs to be square");
	int D = G.size1();
	Vec<Seq<T>> g(D);
	for( int d = 0; d < D; d++ ) {
		for( int e = 0; e < D; e++ ) {
			const Seq<T>& ph = upsample( G(d,e), D );
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
