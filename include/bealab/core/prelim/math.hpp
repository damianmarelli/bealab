/// @file bealab/core/prelim/math.hpp
/// A number of math functions for real and complex numbers.

#ifndef _BEALAB_PRELIM_MATH_
#define	_BEALAB_PRELIM_MATH_

#include <cmath>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/math/special_functions/binomial.hpp>

namespace bealab
{
/// @defgroup prelim_math Math functions
/// A number of math functions, most of them imported from STD.
/// - Classification
///   - isnan
///   - isinf
///   - isfinite
/// - Basic functions
///   - min
///   - max
///   - round
///   - trunc
///   - ceil
///   - floor
///   - mod
/// - Complex functions
///   - abs
///   - arg
///   - real
///   - imag
///   - conj
///   - polar
/// - Trigonometric functions
///   - sin
///   - cos
///   - tan
///   - asin
///   - acos
///   - atan
/// - Hyperbolic functions
///   - sinh
///   - cosh
///   - tanh
///   - asinh
///   - acosh
///   - atanh
/// - Exponential and power functions
///   - exp
///   - log
///   - log2
///   - log10
///   - pow
///   - sqrt
/// - Special functions
///   - sinc
///   - factorial
///   - erf
///   - erfc
///   - erf_inv
///   - erfc_inv
/// @{

/// @name Classification
using std::isnan;
using std::isinf;
using std::isfinite;
/// @}

/// @name Basic functions
using std::min;
using std::max;
#ifdef BEALAB_MACOSX
using ::round;
using ::trunc;
#else
using std::round;
using std::trunc;
#endif
using std::ceil;
using std::floor;

/// x modulo y
inline double mod( double x, double y ) { return x - floor( x / y ) * y; }
/// @}

/// @name Complex functions
using std::abs;
using std::arg;
using std::real;
using std::imag;
using std::conj;
using std::polar;
//inline double abs( const complex& x )  { return std::abs(_complex(x)); }
//inline double arg( const complex& x ) { return std::arg(_complex(x)); }
//inline double real( const complex& x ) { return std::real(_complex(x)); }
//inline double imag( const complex& x ) { return std::imag(_complex(x)); }
//inline complex conj( const complex& x ) { return std::conj(_complex(x)); }
//inline complex polar( double r, double theta=0 ) { return std::polar(r,theta); }
/// @}

/// @name Trigonometric functions
using std::sin;
using std::cos;
using std::tan;
using std::asin;
using std::acos;
using std::atan;
//inline complex sin( const complex& x ) { return std::sin(_complex(x)); }
//inline complex cos( const complex& x ) { return std::cos(_complex(x)); }
//inline complex tan( const complex& x ) { return std::tan(_complex(x)); }
//inline complex asin( const complex& x ) { return std::asin(_complex(x)); }
//inline complex acos( const complex& x ) { return std::acos(_complex(x)); }
//inline complex atan( const complex& x ) { return std::atan(_complex(x)); }
/// @}

/// @name Hyperbolic functions
using std::sinh;
using std::cosh;
using std::tanh;
using std::asinh;
using std::acosh;
using std::atanh;
//inline complex sinh( const complex& x ) { return std::sinh(_complex(x)); }
//inline complex cosh( const complex& x ) { return std::cosh(_complex(x)); }
//inline complex tanh( const complex& x ) { return std::tanh(_complex(x)); }
//inline complex asinh( const complex& x ) { return std::asinh(_complex(x)); }
//inline complex acosh( const complex& x ) { return std::acosh(_complex(x)); }
//inline complex atanh( const complex& x ) { return std::atanh(_complex(x)); }
/// @}

/// @name Exponential and power functions
using std::exp;
using std::log;
#ifdef BEALAB_MACOSX
using ::log2;
#else
using std::log2;
#endif
using std::log10;
//inline complex exp( const complex& x ) { return std::exp(_complex(x)); }
//inline complex log( const complex& x ) { return std::log(_complex(x)); }
//inline complex log2( const complex& x ) { return log(_complex(x))/log(2); }
//inline complex log10( const complex& x ) { return std::log10(_complex(x)); }
using std::pow;
//inline double pow( double x, int y ) { return std::pow( x, y ); }
//inline complex pow( double x, double y )  { return std::pow( _complex(x), y ); }
//inline complex pow( const complex& x, int y ) { return std::pow(_complex(x),y); }
//inline complex pow( const complex& x, double y ) { return std::pow(_complex(x),y); }
//inline complex pow( double x, const complex& y ) { return std::pow(x,_complex(y)); }
using std::sqrt;
//inline complex sqrt( const complex& x ) { return std::sqrt(_complex(x)); }
/// @}

/// @name Special functions

/// Sine cardinal function
//inline double sinc( double x ) { return x == 0 ? 1 : sin( pi * x ) / (pi * x); }
inline double sinc( double x ) { return boost::math::sinc_pi( pi * x ); }

/// Factorial
inline double factorial( int i ) { return boost::math::factorial<double>(i); }

/// Binomial coefficients
inline double binomial_coefficient( int n, int k )
	{ return boost::math::binomial_coefficient<double>(n,k); }

using boost::math::erf;
using boost::math::erfc;
using boost::math::erf_inv;
using boost::math::erfc_inv;

/// @}

/// @}
}
#endif
