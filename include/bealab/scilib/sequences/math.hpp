// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/scilib/sequences/math.hpp
/// Some math operations applied to each entry of the sequence.

#ifndef _BEALAB_SEQUENCES_MATH_
#define	_BEALAB_SEQUENCES_MATH_

namespace bealab
{
/// @defgroup sequences_math Entry-wise functions
/// Some math operations applied to each entry of the sequence.
/// @{

/// Macro to turn a scalar function into a function that applies to each entry.
#define _SEQ_FUNCTION( fun ) \
template<class T> \
auto fun( const sequence<T>& v ) -> sequence<decltype( fun(T()) )> \
{ \
	typedef decltype( fun(T()) ) R; \
	return sequence<R>( fun( v.buffer() ), v.t1() ); \
}

/// @name Roundup functions
_SEQ_FUNCTION( round )
_SEQ_FUNCTION( ceil )
_SEQ_FUNCTION( floor )
_SEQ_FUNCTION( trunc )
_SEQ_FUNCTION( mod )
/// @}

// Complex
_SEQ_FUNCTION( real )
_SEQ_FUNCTION( imag )
_SEQ_FUNCTION( abs )
_SEQ_FUNCTION( arg )
_SEQ_FUNCTION( conj )

// Trigonometric functions
_SEQ_FUNCTION( sin )
_SEQ_FUNCTION( cos )
_SEQ_FUNCTION( tan )
_SEQ_FUNCTION( sinh )
_SEQ_FUNCTION( cosh )
_SEQ_FUNCTION( tanh )

_SEQ_FUNCTION( sqrt )
_SEQ_FUNCTION( exp )
_SEQ_FUNCTION( log )
_SEQ_FUNCTION( log2 )
_SEQ_FUNCTION( log10 )
/// @}

// XXX Faltan las power functions

/// @}
}

#endif
