// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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
