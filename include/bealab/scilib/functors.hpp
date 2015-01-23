// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/scilib/functors.hpp
/// Extends the functionality of std::function

#ifndef _BEALAB_FUNCTORS_
#define	_BEALAB_FUNCTORS_

#include <bealab/core/blas.hpp>
#include <bealab/scilib/sequences.hpp>

namespace bealab
{
/// @defgroup functors Functors
/// Extends the functionality of std::function
/// @{

/// Sum of functors
template<class R, class... A>
function<R(A)...> operator+( const function<R(A...)>& f1, const function<R(A...)>& f2 )
{
	return [f1,f2]( A... a ) { return f1(a...) + f2(a...); };
}

/// Multiplication of functors
template<class R, class... A>
function<R(A)...> operator*( const function<R(A...)>& f1, const function<R(A...)>& f2 )
{
	return [f1,f2]( A... a ) { return f1(a...) * f2(a...); };
}

/// Pre-multiplication by a scalar
template<class R, class... A, class S>
function<R(A)...> operator*( const S& c, const function<R(A...)>& f1 )
{
	return [f1,c]( A... a ) { return c * f1(a...); };
}

/// Post-multiplication by a scalar
template<class R, class... A, class S>
function<R(A)...> operator*( const function<R(A...)>& f1, const S& c )
{
	return c * f1;
}

/// Inverse of a functor
function<double(double)> inv( const function<double(double)>& fun, double lo, double hi );

/// @}
}
#endif
