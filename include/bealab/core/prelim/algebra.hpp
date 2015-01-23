// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/core/prelim/algebra.hpp
/// Operation for all basic types used to build algebras.

#ifndef _BEALAB_PRELIM_ALGEBRA_
#define	_BEALAB_PRELIM_ALGEBRA_

namespace bealab
{
/// @defgroup prelim_algebra Abstract algebra
/// Operations to build algebraic structures using double and complex types.
/// @{

/// @name Real numbers (double)

/// Adjoint operation
inline
double adjoint( double x ) { return x; }

/// Inverse
inline
double inv( double x ) { return 1/x; }

/// Norm
inline
double norm( double x, int p=2 ) { return p!=0 ? abs(x) : (x==0 ? 0 : 1); }

/// Inner product
inline
double inner_prod( double x, double y ) { return x*y; }

/// Outer product
inline
double outer_prod( double x, double y ) { return x*y; }
/// @}

/// @name Complex numbers

/// Adjoint operation
inline
complex adjoint( const complex &x ) { return conj(x); }

/// Inverse
inline
complex inv( const complex &x ) { return 1/x; }

/// Norm
inline
double norm( const complex &x, int p=2 ) { return p!=0 ? abs(x) : (abs(x)==0 ? 0 : 1); }

/// Inner product
inline
complex inner_prod( const complex &x, const complex &y ) { return x*conj(y); }

/// Outer product
inline
complex outer_prod( const complex &x, const complex &y ) { return x*conj(y); }
/// @}

/// @}
}
#endif
