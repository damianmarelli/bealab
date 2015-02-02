/// @file bealab/core/prelim/complex.hpp
/// Preliminary file.

#ifndef _BEALAB_PRELIM_COMPLEX_
#define	_BEALAB_PRELIM_COMPLEX_

#include <complex>

namespace bealab
{
/// @defgroup prelim_complex Complex numbers
/// Define complex as std::complex<double>, and add support for algebraic
/// operations with integers and booleans.
/// @{

/// Define complex as std::complex<double>
typedef std::complex<double> complex;

/*
/// Class of complex numbers.
/// Complex class supporting automatic casting to double and algebraic operations with integers.
class complex : public _complex {
public:

	/// Inherit the constructors from base
	INHERIT_CONSTRUCTORS( complex, _complex );

	/// Inherit the assignments from base
	INHERIT_ASSIGNMENTS( complex, _complex );

//	operator double() const
//	{
//		return this->real();
//	}																			///< Cast operator
};
*/

/// Display
template<class A, class B>
std::basic_ostream<A,B> &operator<<( std::basic_ostream<A,B> &os, const complex &x )
{
	os << x.real();
	if( x.imag() >= 0 )
		os << "+";
	os << x.imag() << "i";
    return os;
}

/// @name Arithmetic operations with bool
inline complex operator+( bool a, const complex& z ) { return (double)a+z; }
inline complex operator+( const complex& z, bool a ) { return z+(double)a; }
inline complex operator-( bool a, const complex& z ) { return (double)a-z; }
inline complex operator-( const complex& z, bool a ) { return z-(double)a; }
inline complex operator*( bool a, const complex& z ) { return (double)a*z; }
inline complex operator*( const complex& z, bool a ) { return z*(double)a; }
inline complex operator/( bool a, const complex& z ) { return (double)a/z; }
inline complex operator/( const complex &z, bool a ) { return z/(double)a; }
/// @}

/// @name Arithmetic operations with int
inline complex operator+( int a, const complex& z ) { return (double)a+z; }
inline complex operator+( const complex& z, int a ) { return z+(double)a; }
inline complex operator-( int a, const complex& z ) { return (double)a-z; }
inline complex operator-( const complex& z, int a ) { return z-(double)a; }
inline complex operator*( int a, const complex& z ) { return (double)a*z; }
inline complex operator*( const complex& z, int a ) { return z*(double)a; }
inline complex operator/( int a, const complex& z ) { return (double)a/z; }
inline complex operator/( const complex &z, int a ) { return z/(double)a; }
/// @}

/// @}
}
#endif
