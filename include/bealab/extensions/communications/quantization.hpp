// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/extensions/communications/quantization.hpp
/// Routines for analysis and implementation of scalar and vector quantizers.

#ifndef _BEALAB_QUANTIZATION_
#define	_BEALAB_QUANTIZATION_

#include <bealab/core/blas.hpp>
#include <bealab/scilib/sequences.hpp>
#include <bealab/extensions/control/linsys.hpp>

namespace bealab
{
namespace comms
{
//------------------------------------------------------------------------------
/// @defgroup quantization Quantization
/// Routines for analysis and implementation of scalar and vector quantizers.
/// @{

class quantizer{
public:
	rvec bnd;
	rvec dict;
	quantizer(){};
	quantizer( const rvec& );
	int encode( double ) const;
	iseq encode( const rseq& ) const;
	double decode( int ) const;
	rseq decode( const iseq& ) const;
	double operator()( double ) const;
	rseq operator()( const rseq& ) const;
};

class vector_quantizer{
public:
	Vec<rvec> dict;
	vector_quantizer(){};
	vector_quantizer( const Vec<rvec> &dict_ ) : dict(dict_) {};
	vector_quantizer( const Vec<rvec>&, int );
	void train( const Vec<rvec>&, int );
	int encode( const rvec& ) const;
	iseq encode( const Vec<rseq>& ) const;
	iseq encode( const Seq<rvec>& ) const;
	rvec decode( int ) const;
	Seq<rvec> decode( const iseq& ) const;
	rvec operator()( const rvec& ) const;
	Vec<rseq> operator()( const Vec<rseq> &x ) const;
	Seq<rvec> operator()( const Seq<rvec> &x ) const;
};

//------------------------------------------------------------------------------
/// @defgroup analysis Analysis
/// Functions for quantization analysis.
/// @{

/// @name For a generic distribution
double I0( const function<double(double)>&, double, double );
double I1( const function<double(double)>&, double, double );
double I2( const function<double(double)>&, double, double );
double qerror( const function<double(double)>&, const quantizer& );
double qpower( const function<double(double)>&, const quantizer& );
double cmean( const function<double(double)>&, double, double );
/// @}

/// @name For the normal distribution
double I0( double, double, double, double );
double I1( double, double, double, double );
double I2( double, double, double, double );
double qerror( double, double, const quantizer& );
double qpower( double, double, const quantizer& );
double cmean( double, double, double, double );
/// @}

/// @name Rate distortion functions
double distortion_rate_function( double R, const function<double(double)>& phi );
/// @}

/// @}

//------------------------------------------------------------------------------
/// @defgroup implementation Implementation
/// Functions for quantizing.
/// @{

/// @name Quantizer design
rvec qboundaries( const rvec& );
rvec qdict( const function<double(double)>&, const rvec& );
quantizer lloyd( const function<double(double)>&, const rvec& );
quantizer lloyd( double, double, int );
rvec waterfilling( const rvec&, int );
Vec<rseq> KL_filterbank( const rseq &, int, rvec* =NULL, double=1e-4 );
/// @}

// LPC scalar
//rarma lpc_predictor( const rseq&, int );
//quantizer lpc_quantizer( int, const rseq&, const rseq&, int );
//rseq lpc( const rseq&, rarma, const quantizer& );

/// @name Linear predictive quantization
control::arma<rmat,rvec> lpc_predictor( const Seq<rmat> &Fx_, int N );
vector_quantizer lpc_quantizer_openloop( int L, const function<Seq<rvec>(Seq<rvec>)>& P,
		const Seq<rmat> &Fx, int K );
vector_quantizer lpc_quantizer_closedloop( int L, const function<Seq<rvec>(Seq<rvec>)>& P,
		const Seq<rmat> &Fx, int K );
Seq<rvec> lpc( const Seq<rvec> &x, const function<rvec(rvec)>& p, const vector_quantizer &q );
Vec<rseq> lpc( const Vec<rseq> &x, const Mat<rseq> &Fx, int P, int L, int K,
		const string& ="closed loop");
/// @}

/// @}

/// @}

}
}
#endif
