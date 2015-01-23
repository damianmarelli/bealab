// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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

}
}
#endif
