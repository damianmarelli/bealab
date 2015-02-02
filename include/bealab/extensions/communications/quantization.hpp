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
	vec<rvec> dict;
	vector_quantizer(){};
	vector_quantizer( const vec<rvec> &dict_ ) : dict(dict_) {};
	vector_quantizer( const vec<rvec>&, int );
	void train( const vec<rvec>&, int );
	int encode( const rvec& ) const;
	iseq encode( const vec<rseq>& ) const;
	iseq encode( const sequence<rvec>& ) const;
	rvec decode( int ) const;
	sequence<rvec> decode( const iseq& ) const;
	rvec operator()( const rvec& ) const;
	vec<rseq> operator()( const vec<rseq> &x ) const;
	sequence<rvec> operator()( const sequence<rvec> &x ) const;
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
vec<rseq> KL_filterbank( const rseq &, int, rvec* =NULL, double=1e-4 );
/// @}

// LPC scalar
//rarma lpc_predictor( const rseq&, int );
//quantizer lpc_quantizer( int, const rseq&, const rseq&, int );
//rseq lpc( const rseq&, rarma, const quantizer& );

/// @name Linear predictive quantization
control::arma<rmat,rvec> lpc_predictor( const sequence<rmat> &Fx_, int N );
vector_quantizer lpc_quantizer_openloop( int L, const function<sequence<rvec>(sequence<rvec>)>& P,
		const sequence<rmat> &Fx, int K );
vector_quantizer lpc_quantizer_closedloop( int L, const function<sequence<rvec>(sequence<rvec>)>& P,
		const sequence<rmat> &Fx, int K );
sequence<rvec> lpc( const sequence<rvec> &x, const function<rvec(rvec)>& p, const vector_quantizer &q );
vec<rseq> lpc( const vec<rseq> &x, const mat<rseq> &Fx, int P, int L, int K,
		const string& ="closed loop");
/// @}

/// @}

/// @}

}
}
#endif
