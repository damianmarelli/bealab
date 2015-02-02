/// @file bealab/core/blas/entrywise.hpp
/// Functions applied to each entry of a vector/matrix.

#ifndef _BEALAB_BLAS_ENTRYWISE_
#define	_BEALAB_BLAS_ENTRYWISE_

#include <bealab/core/blas/matrix.hpp>

namespace bealab
{
/// @defgroup blas_entrywise Entry-wise functions
/// Functions applied to each entry of a vector/matrix.
/// @{

/// Macro to turn a scalar functions into a vector one and a matrix one.
/// This macro converts a scalar function into two functions. One applies the
/// scalar function to each entry of a vector, and the other does the same
/// with a matrix.
#define _ENTRYWISE_FUNCTION( fun ) \
template<class E> \
vec<decltype( fun(typename E::value_type()) )> fun( const vector_interface<E>& v ) \
{ \
	typedef decltype( fun(typename E::value_type()) ) R; \
	int I = v.size(); \
	vec<R> r(I); \
	for( int i = 0;  i < I; i++ ) \
		r(i) = fun(v(i)); \
	return r; \
} \
\
template<class E> \
mat<decltype( fun(typename E::value_type()) )> fun( const matrix_interface<E>& A ) \
{ \
	typedef decltype( fun(typename E::value_type()) ) R; \
	int I = A.size1(); \
	int J = A.size2(); \
	mat<R> B(I,J); \
	for( int i = 0;  i < I; i++ ) \
		for( int j = 0;  j < J; j++ ) \
			B(i,j) = fun(A(i,j)); \
	return B; \
}

/// @name Entry-wise basic functions
_ENTRYWISE_FUNCTION( round );
_ENTRYWISE_FUNCTION( ceil );
_ENTRYWISE_FUNCTION( floor );
_ENTRYWISE_FUNCTION( trunc );

template<class T, class S>
vec<decltype(mod(typename T::value_type(),S()))>
mod( const vector_interface<T>& x, S a )
{
	typedef decltype(mod(typename T::value_type(),S())) R;
	int N = x.size();
	vec<R> rv(N);
	for( int n = 0; n < N; n++ )
		rv(n) = mod( x(n), a );
	return rv;
}

template<class T, class S>
mat<decltype(mod(typename T::value_type(),S()))>
mod( const matrix_interface<T>& x, S a )
{
	typedef decltype(mod(typename T::value_type(),S())) R;
	int I = x.size1();
	int J = x.size2();
	mat<R> rv(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			rv(i,j) = mod( x(i,j), a );
	return rv;
}
/// @}

/// @name Entry-wise complex funtions
_ENTRYWISE_FUNCTION( real );
_ENTRYWISE_FUNCTION( imag );
_ENTRYWISE_FUNCTION( abs );
_ENTRYWISE_FUNCTION( arg );
_ENTRYWISE_FUNCTION( conj );
/// @}

/// @name Entry-wise trigonometric functions
_ENTRYWISE_FUNCTION( sin );
_ENTRYWISE_FUNCTION( cos );
_ENTRYWISE_FUNCTION( tan );
_ENTRYWISE_FUNCTION( asin );
_ENTRYWISE_FUNCTION( acos );
_ENTRYWISE_FUNCTION( atan );
/// @}

/// @name Entry-wise hyperbolic functions
_ENTRYWISE_FUNCTION( sinh );
_ENTRYWISE_FUNCTION( cosh );
_ENTRYWISE_FUNCTION( tanh );
_ENTRYWISE_FUNCTION( asinh );
_ENTRYWISE_FUNCTION( acosh );
_ENTRYWISE_FUNCTION( atanh );
/// @}

/// @name Entry-wise exponential functions
_ENTRYWISE_FUNCTION( exp );
_ENTRYWISE_FUNCTION( log );
_ENTRYWISE_FUNCTION( log2 );
_ENTRYWISE_FUNCTION( log10 );
/// @}

/// @name Entry-wise power functions
_ENTRYWISE_FUNCTION( sqrt );

template<class T, class S>
vec<decltype(pow(typename T::value_type(),typename S::value_type()))>
pow( const vector_interface<T>& x, const vector_interface<S>& y )
{
	assert( x.size() == y.size() );
	typedef decltype(pow(typename T::value_type(),typename S::value_type())) R;
	int N = x.size();
	vec<R> rv(N);
	for( int n = 0; n < N; n++ )
		rv(n) = pow( x(n), y(n) );
	return rv;
}

template<class T>
vec<decltype(pow(typename T::value_type(),double(0)))>
pow( const vector_interface<T>& x, double a )
{
	return pow( x, a*ones(x.size()) );
}

template<class T>
vec<decltype(pow(typename T::value_type(),complex(0,0)))>
pow( const vector_interface<T>& x, complex a )
{
	return pow( x, a*ones(x.size()) );
}

template<class T>
vec<decltype(pow(double(0),typename T::value_type()))>
pow( double a, const vector_interface<T>& x )
{
	return pow( a*ones(x.size()), x );
}

template<class T>
vec<decltype(pow(complex(0,0),typename T::value_type()))>
pow( complex a, const vector_interface<T>& x )
{
	return pow( a*ones(x.size()), x );
}

template<class T, class S>
mat<decltype(pow(typename T::value_type(),typename S::value_type()))>
pow( const matrix_interface<T>& x, const matrix_interface<S>& y )
{
	assert( x.size1() == y.size1() );
	assert( x.size2() == y.size2() );
	typedef decltype(pow(typename T::value_type(),typename S::value_type())) R;
	int I = x.size1();
	int J = x.size2();
	mat<R> rv(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			rv(i,j) = pow( x(i,j), y(i,j) );
	return rv;
}

template<class T>
mat<decltype(pow(typename T::value_type(),double(0)))>
pow( const matrix_interface<T>& x, double a )
{
	return pow( x, a*ones( x.size1(),
							  x.size2() ) );
}

template<class T>
mat<decltype(pow(typename T::value_type(),complex(0,0)))>
pow( const matrix_interface<T>& x, complex a )
{
	return pow( x, a*ones( x.size1(),
							  x.size2() ) );
}

template<class T>
mat<decltype(pow(double(0),typename T::value_type()))>
pow( double a, const matrix_interface<T>& x )
{
	return pow( a*ones( x.size1(),
						x.size2() ), x );
}

template<class T>
mat<decltype(pow(complex(0,0),typename T::value_type()))>
pow( complex a, const matrix_interface<T>& x )
{
	return pow( a*ones( x.size1(),
						x.size2() ), x );
}
/// @}

/// Entry-wise function made from a scalar function.
template<class F>
class entrywise_function {

	F fun;																		///< Scalar function to be applied

public:

	/// Constructor
	entrywise_function( const F &f ) : fun(f) {}

	/// Apply the scalar function entry-wise
	template<class T>
	auto operator()( const T &x ) -> decltype(entrywise(fun,x))
	{
		return entrywise( fun, x );
	}
};

/// Conversion a scalar function into an entry-wise function.
/// For each type for which you want to define entry-wise functions, you need to
/// define the function template<class FUN> type entrywise( FUN, type )
template<class F>
entrywise_function<F> entrywise( F fun )
{
	return fun;
}

/// @}
}
#endif
