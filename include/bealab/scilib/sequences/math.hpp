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
auto fun( const Seq<T>& v ) -> Seq<decltype( fun(T()) )> \
{ \
	typedef decltype( fun(T()) ) R; \
	return Seq<R>( fun( v.vec() ), v.t1() ); \
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
