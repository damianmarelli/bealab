// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/extensions/signal/tfanalysis.hpp
/// Routines for time-frequency analysis.

#ifndef _BEALAB_EXTENSIONS_TFANALYSIS_
#define	_BEALAB_EXTENSIONS_TFANALYSIS_

#include <bealab/extensions/signal/mrsigpro.hpp>

namespace bealab
{
namespace signal
{
//------------------------------------------------------------------------------
/// @defgroup tfanalysis Time-frequency analysis
/// Routines for time-frequency analysis.
/// @{

/// @name Analysis and synthesis

/// Analysis using a generic filterbank
template <class T, class S>
Vec<Seq<decltype(S()*T())>> analysis( const Seq<T> &x, const Vec<Seq<S>> &h, int D )
{
	const Vec<Seq<T>>& Xpp = fb2pp_signal( x, D );
	const Mat<Seq<S>>& H   = fb2pp_filterbank( h, D );
	return H * Xpp;
}

/// Synthesis using a generic filterbank
template <class T, class S>
Seq<decltype(S()*T())> synthesis( const Vec<Seq<T>> &X, const Vec<Seq<S>> &h, int D )
{
	typedef decltype(S()*T()) R;
	const Mat<Seq<S>>& H   = fb2pp_filterbank( h, D );
	const Vec<Seq<R>>& Xpp = H.A() * X;
	return pp2fb_signal( Xpp );
}

/// Analysis using a DFT filterbank
template<class T, class S>
Vec<cseq> analysis ( const Seq<T>& x, const Seq<S>& h0, int M, int D )
{
	// If the input is empty
	if( x.size() == 0 )
		return Vec<cseq>(M);

	// Constants
	int h1 = h0.t1();
	int h2 = h0.t2();
	int m1 = floor( h1     / (double)M );
	int m2 = floor( h2     / (double)M );
	int d1 = floor( h1     / (double)D );
	int d2 = floor( h2     / (double)D );
	int k1 = ceil( x.t1()  / (double)D ) + d1;
	int k2 = ceil( x.t2()  / (double)D ) + d2;
	int K  = k2-k1+1;

	// Type of the product
	typedef decltype(T()*S()) R;

	// Prepare output signals
	Vec<cseq> X(M);
	for( int m = 0; m < M; m++ )
		X(m) = cseq( zeros(K), k1 );

	// For each subband sample
	for( int k = k1; k <= k2; k++ ) {

		// Take a flipped piece of the input
		const Vec<T>& u1 = x.trunc( k*D-h2, k*D-h1 ).vec();
		const Seq<T>& u2 = Seq<T>( flip(u1), h1 );

		// Multiply it with the window
		const Seq<R>& u3 = element_prod( u2, h0 );

		// Overlap-add u3
		Vec<R> u4 = zeros(M);
		for( int m = m1; m <= m2; m++ )
			u4 += u3.trunc( m*M, (m+1)*M-1 ).vec();

		// Take the inverse Fourier transform
		const cvec& u5 = M*idft( u4 );

		// Assign it to the output signals
		for( int m = 0; m < M; m++ )
			X(m)(k) = u5(m);
	}

	return X;
}

/// Synthesis using a DFT filterbank
template<class T, class S>
cseq synthesis ( const Vec<Seq<T>>& X, const Seq<S>& h0, int D )
{
	// Constants
	int M  = X.size();
	int h1 = h0.t1();
	int h2 = h0.t2();
	int m1 = floor( h1     / (double)M );
	int m2 = floor( h2     / (double)M );
	int k1 = min( entrywise([](const Seq<T>&x) {return x.t1();})(X) );
	int k2 = max( entrywise([](const Seq<T>&x) {return x.t2();})(X) );

	// Output signal
	int y1 = k1*D-h2;
	int y2 = k2*D-h1;
	int Ly = y2 - y1 + 1;
	cseq y( Ly, y1 );

	// For each subband sample
	for( int k = k1; k <= k2; k++ ) {

		// Take the subband sammple in a vector
		Vec<T> u1(M);
		for( int m = 0; m < M; m++ )
			if( k >= X(m).t1() && k <= X(m).t2() )
				u1(m) = X(m)(k);
			else
				u1(m) = 0;

		// Take the Fourier transform
		const cvec& u2 = dft( u1 );

		// For a vector of m2-m1+1 repeated copies
		cvec u3((m2-m1+1)*M);
		for( int m = m1; m <= m2; m++ )
			u3( range(m*M,(m+1)*M) ) = u2;
		const cseq& u4 = cseq( u3, m1*M );

		// Multiply it with the window
		const cseq& u5 = element_prod( u4, h0 );

		// Overlap-add it to the output
		const cseq& u6 = u5.timereverse();
		int t1 = k*D - h2;
		int t2 = k*D - h1;
		y.vec()( range( t1-y1, t2-y1+1 ) ) += u6.vec();
	}

	return y;
}
/// @}

/// @name Filterbank design

/// Build a DFT filterbank with the window h0
Vec<cseq> filterbank_dtft( const cseq& h0, int M );

/// Obtain the best possible dual of the window h, with support in [t1,t2]
cseq quasidual_window( const cseq& h, int M, int D, int t1, int t2 );

/// Obtain the canonical dual of the window h, within tolerance tol on the reconstruction error
cseq connonical_dual_window( const cseq& h0, int M, int D, double tol );

/// Compute the reconstruction error of a pair window/dual
double dual_window_error( const cseq& h0, const cseq& f0, int M, int D );

/// @}

/// @name Other

/// Builld a spectrogram of the signal x, using M frequencies, downsampling
/// factor D and a Hamming window of length N
cmat spectrogram( rseq x, int M, int D, int N );
/// @}

/// @}
}
}

#endif
