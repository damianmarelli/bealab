// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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
