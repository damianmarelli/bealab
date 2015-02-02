/// @file bealab/core/blas/reductions.hpp
/// Functionals applied to a vector/matrix, or to each row/column of a matrix.

#ifndef _BEALAB_BLAS_REDUCTIONS_
#define	_BEALAB_BLAS_REDUCTIONS_

#include <bealab/core/blas/matrix.hpp>

namespace bealab
{
/// @defgroup blas_reductions Reductions
/// Functionals applied to a vector/matrix, or to each row/column of a matrix.
/// @{

/// @name Macro to convert a functional into a reduction

/// Macro to turn a vector functional into a row/column-wise reduction of a matrix.
#define _ROWCOLWISE_FUNCTIONAL( fun ) \
template<class T> \
vec<typename T::value_type> fun##_rowwise( const matrix_interface<T>& A ) \
{ \
	typedef decltype( fun( A.row(0) ) ) R; \
	int N = A.size1(); \
	vec<R> rv(N); \
	for( int n = 0; n < N; n++ ) \
		rv(n) = fun( A.row(n) ); \
	return rv; \
} \
\
template<class T> \
vec<typename T::value_type> fun##_columnwise( const matrix_interface<T>& A ) \
{ \
	typedef decltype( fun( A.column(0) ) ) R; \
	int N = A.size2(); \
	vec<R> rv(N); \
	for( int n = 0; n < N; n++ ) \
		rv(n) = fun( A.column(n) ); \
	return rv; \
}

/// Macro to turn a vector functional into a matrix functional.
#define _MATRIX_FUNCTIONAL( fun ) \
template<class T> \
typename T::value_type fun( const matrix_interface<T>& A ) \
{ \
	return fun( fun##_columnwise( A ) ); \
}

/// @}

/// @name Find an entry
inline
ivec find( const bvec &x )
{
	int N = x.size();
	int K = sum(ivec(x));
	ivec rv(K);
	for( int n = 0, k = 0; n < N; n++ )
		if( x(n) == 1 )
			rv(k++) = n;
	return rv;
}

template<class T>
int max_index( const vector_interface<T>& x )
{
	assert( x.size() > 0 );
	typedef typename T::value_type R;
	R mv = x(0);
	int mi = 0;
	int N = x.size();
	for( int n = 1; n < N; n++ )
		if( mv < x(n) ) {
			mv = x(n);
			mi = n;
		}
	return mi;
}

template<class T>
int min_index( const vector_interface<T>& x )
{
	return max_index(-x);
}

_ROWCOLWISE_FUNCTIONAL( find )
_ROWCOLWISE_FUNCTIONAL( max_index )
_ROWCOLWISE_FUNCTIONAL( min_index )
/// @}

/// @name Max/min
template<class T>
typename T::value_type max( const vector_interface<T>& x )
{
	if( x.size() == 0 )
		return -std::numeric_limits<typename T::value_type>::infinity();
	return x(max_index(x));
}

template<class T>
typename T::value_type min( const vector_interface<T>& x )
{
	return -max(-x);
}

_ROWCOLWISE_FUNCTIONAL( max )
_ROWCOLWISE_FUNCTIONAL( min )
_MATRIX_FUNCTIONAL( max )
_MATRIX_FUNCTIONAL( min )
/// @}

/// @name Sum, difference and product of entries
template<class T>
typename T::value_type sum( const vector_interface<T>& x )
{
	typedef typename T::value_type R;
	if( x.size() == 0 )
		return R();
	R rv = x(0);
	int N = x.size();
	for( int n = 1; n < N; n++ )
		rv += x(n);
	return rv;
}

template<class T>
vec<typename T::value_type> diff( const vector_interface<T>& x )
{
	typedef typename T::value_type R;
	if( x.size() <= 1 )
		return vec<R>();
	int N = x.size();
	vec<R> rv(N-1);
	for( int n = 1; n < N; n++ )
		rv(n-1) = x(n) - x(n-1);
	return rv;
}

template<class T>
typename T::value_type prod( const vector_interface<T>& x )
{
	assert( x.size() > 0 );
	typedef typename T::value_type R;
	R rv = x(0);
	int N = x.size();
	for( int n = 1; n < N; n++ )
		rv *= x(n);
	return rv;
}

_ROWCOLWISE_FUNCTIONAL( sum )
_ROWCOLWISE_FUNCTIONAL( prod )
_MATRIX_FUNCTIONAL( sum )
_MATRIX_FUNCTIONAL( prod )
/// @}

/// @}
}
#endif
