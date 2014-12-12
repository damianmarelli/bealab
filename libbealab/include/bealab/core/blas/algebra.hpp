/// @file bealab/core/blas/algebra.hpp
/// Algebraic operations

#ifndef _BEALAB_BLAS_ALGEBRA_
#define	_BEALAB_BLAS_ALGEBRA_

#include <bealab/core/blas/reductions.hpp>
#include <bealab/core/blas/entrywise.hpp>

namespace bealab
{
/// @defgroup blas_algebra Algebra
/// Algebraic operations
/// @{

/*
 * Do this: If double or complex, use BLAS. Otherwise, use brute force.
 */

/// Check for scalar types (i.e., bool, int, double complex)
//template< class T >
//struct is_scalar : std::integral_constant<
//	   bool,
//	   std::is_same<bool, typename std::remove_cv<T>::type>::value  ||
//	   std::is_same<int, typename std::remove_cv<T>::type>::value  ||
//	   std::is_same<double, typename std::remove_cv<T>::type>::value  ||
//	   std::is_same<complex, typename std::remove_cv<T>::type>::value
//   > {};

/// @name Vector operations

/// vector + vector
template< class E1, class E2,
		class R = decltype(typename E1::value_type()+typename E2::value_type()) >
Vec<R> operator+( const vector_interface<E1>& x, const vector_interface<E2>& y )
{
	assert( x.size() == y.size() );
	int I = x.size();
	Vec<R> z(I);
	for( int i = 0; i < I; i++ )
		z(i) = x(i) + y(i);
	return z;
}
//template< class E1, class E2,
//	class R = decltype( declval<E1>() + declval<E2>() ) >
//vector_interface<R> operator+( const vector_interface<E1>& x, const vector_interface<E2>& y )
//{
//	return dynamic_cast<const E1&>(x) + dynamic_cast<const E2&>(y);
//}

/// -vector
template<class E, class R = typename E::value_type>
Vec<R> operator-( const vector_interface<E>& x )
{
	int I = x.size();
	Vec<R> z(I);
	for( int i = 0; i < I; i++ )
		z(i) = -x(i);
	return z;
}
//template<class E,
//	class R = decltype( -declval<E>() ) >
//vector_interface<R> operator-( const vector_interface<E>& x )
//{
//	return -dynamic_cast<const E&>(x);
//}

/// vector - vector
template<class E1, class E2,
	class R = decltype(typename E1::value_type()-typename E2::value_type())>
Vec<R> operator-( const vector_interface<E1>& x, const vector_interface<E2>& y )
{
	assert( x.size() == y.size() );
	int I = x.size();
	Vec<R> z(I);
	for( int i = 0; i < I; i++ )
		z(i) = x(i) - y(i);
	return z;
}
//template<class E1, class E2,
//	class R = decltype( declval<E1>() - declval<E2>() ) >
//vector_interface<R> operator-( const vector_interface<E1>& x, const vector_interface<E2>& y )
//{
//	return dynamic_cast<const E1&>(x) - dynamic_cast<const E2&>(y);
//}

/// vector * scalar
template< class E, class T,
	class R = decltype( typename E::value_type()*T() ) >
Vec<R> operator*( const vector_interface<E>& x, const T& y )
{
	int I = x.size();
	Vec<R> z(I);
	for( int i = 0; i < I; i++ )
		z(i) = x(i) * y;
	return z;
}
//template<class E, class T,
//	class R = decltype( declval<E>() * declval<T>() ) >
//vector_interface<R> operator*( const vector_interface<E> &x, const T &y )
//{
//	return dynamic_cast<const E&>(x) * y;
//}

/// scalar * vector
template< class T, class E,
	class R = decltype( T()*typename E::value_type() )>
Vec<R> operator*( const T& x, const vector_interface<E>& y )
{
	int I = y.size();
	Vec<R> z(I);
	for( int i = 0; i < I; i++ )
		z(i) = x * y(i);
	return z;
}
//template<class T, class E,
//	class R = decltype( declval<T>() * declval<E>() ) >
//vector_interface<R> operator*( const T& x, const vector_interface<E> &y )
//{
//	return x * dynamic_cast<const E&>(y);
//}

/// vector / scalar
template< class E, class T,
	class R = decltype( typename E::value_type()/T() ) >
Vec<R> operator/( const vector_interface<E>& x, const T& y )
{
	int I = x.size();
	Vec<R> z(I);
	for( int i = 0; i < I; i++ )
		z(i) = x(i) / y;
	return z;
}
/// @}

//------------------------------------------------------------------------------
/// @name Matrix operations

/// matrix + matrix
template< class E1, class E2,
	class R = decltype( typename E1::value_type() + typename E2::value_type() ) >
Mat<R> operator+( const matrix_interface<E1>& x, const matrix_interface<E2>& y )
{
	assert( x.size1() == y.size1() );
	assert( x.size2() == y.size2() );
	int I = x.size1();
	int J = x.size2();
	Mat<R> z(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			z(i,j) = x(i,j) + y(i,j);
	return z;
}
//template< class E1, class E2,
//	class R = decltype( declval<E1>() + declval<E2>() ) >
//matrix_interface<R> operator+( const matrix_interface<E1>& x, const matrix_interface<E2>& y )
//{
//	return dynamic_cast<const E1&>(x) + dynamic_cast<const E2&>(y);
//}

/// -matrix
template< class E, class R = typename E::value_type >
Mat<R> operator-( const matrix_interface<E>& x )
{
	int I = x.size1();
	int J = x.size2();
	Mat<R> z(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			z(i,j) = -x(i,j);
	return z;
}
//template< class E,
//	class R = decltype( -declval<E>() ) >
//matrix_interface<R> operator-( const matrix_interface<E>& x )
//{
//	return -dynamic_cast<const E&>(x);
//}

/// matrix - matrix
template< class E1, class E2,
	class R = decltype( typename E1::value_type() - typename E2::value_type() ) >
Mat<R> operator-( const matrix_interface<E1>& x, const matrix_interface<E2>& y )
{
	assert( x.size1() == y.size1() );
	assert( x.size2() == y.size2() );
	int I = x.size1();
	int J = x.size2();
	Mat<R> z(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			z(i,j) = x(i,j) - y(i,j);
	return z;
}
//template< class E1, class E2,
//	class R = decltype( declval<E1>() - declval<E2>() ) >
//matrix_interface<R> operator-( const matrix_interface<E1>& x, const matrix_interface<E2>& y )
//{
//	return dynamic_cast<const E1&>(x) - dynamic_cast<const E2&>(y);
//}

/// matrix * scalar
template< class E, class T,
	class R = decltype( typename E::value_type()*T() ) >
Mat<R> operator*( const matrix_interface<E>& x, const T& y )
{
	int I = x.size1();
	int J = x.size2();
	Mat<R> z(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			z(i,j) = x(i,j) * y;
	return z;
}

/// scalar * matrix
template< class T, class E,
	class R = decltype( T() * typename E::value_type() ) >
Mat<R> operator*( const T& x, const matrix_interface<E>& y )
{
	int I = y.size1();
	int J = y.size2();
	Mat<R> z(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			z(i,j) = x * y(i,j);
	return z;
}

/// matrix / scalar
template< class E, class T,
	class R = decltype( typename E::value_type() / T() ) >
Mat<R> operator/( const matrix_interface<E>& x, const T& y )
{
	int I = x.size1();
	int J = x.size2();
	Mat<R> z(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			z(i,j) = x(i,j) / y;
	return z;
}
/// @}

//------------------------------------------------------------------------------
/// @name Vector-matrix products

/// matrix * matrix
template< class E1, class E2,
	class R = decltype( typename E1::value_type() * typename E2::value_type() ) >
Mat<R> operator*( const matrix_interface<E1>& A, const matrix_interface<E2>& B )
{
//	return ublas::prod(A,B);
	assert( A.size2() == B.size1() );
	int I = A.size1();
	int J = B.size2();
	int K = A.size2();
	Mat<R> M(I,J);
	for( int i = 0; i < I; i++ ) {
		for( int j = 0; j < J; j++ ) {
			M(i,j) = R();
			for( int k = 0; k < K; k++ )
				if( k == 0 )
					M(i,j)  = A(i,k) * B(k,j);
				else
					M(i,j) += A(i,k) * B(k,j);
		}
	}
	return M;
}

/// matrix * vector
template< class E1, class E2,
	class R = decltype( typename E1::value_type() * typename E2::value_type() ) >
Vec<R> operator*( const matrix_interface<E1>& A, const vector_interface<E2>& x )
{
	assert( A.size2() == x.size() );
	int I = A.size1();
	int K = A.size2();
	Vec<R> v(I);
	for( int i = 0; i < I; i++ ) {
		v(i) = R();
		for( int k = 0; k < K; k++ )
			if( k == 0 )
				v(i)  = A(i,k) * x(k);
			else
				v(i) += A(i,k) * x(k);
	}
	return v;
}

/// vector * matrix
template< class E1, class E2,
	class R = decltype( typename E1::value_type() * typename E2::value_type() ) >
Vec<R> operator*( const vector_interface<E1>& x, const matrix_interface<E2>& B )
{
	assert( x.size() == B.size1() );
	int J = B.size2();
	int K = x.size();
	Vec<R> v(J);
	for( int j = 0; j < J; j++ ) {
		v(j) = R();
		for( int k = 1; k < K; k++ )
			if( k == 0 )
				v(j)  = x(k) * B(k,j);
			else
				v(j) += x(k) * B(k,j);
	}
	return v;
}
/// @}

//------------------------------------------------------------------------------

/// @name Element-wise product/division

/// Vector-vector element-wise product
template<class T, class S,
	class R = decltype(typename T::value_type()*typename S::value_type())>
Vec<R> element_prod( const vector_interface<T> &x, const vector_interface<S> &y )
{
	assert( x.size() == y.size() );
	int I = x.size();
	Vec<R> r(I);
	for( int i = 0; i < I; i++ )
		r(i) = x(i) * y(i);
	return r;
}

/// Vector element-wise inversion
template<class T, class R = typename T::value_type>
Vec<R> element_inv( const vector_interface<T> &x )
{
	int I = x.size();
	Vec<R> r(I);
	for( int i = 0; i < I; i++ )
		r(i) = inv( x(i) );
	return r;
}

/// Vector-vector element-wise division
template<class T, class S,
	class R = decltype(typename T::value_type()*typename S::value_type())>
Vec<R> element_div( const vector_interface<T> &x, const vector_interface<S> &y )
{
	assert( x.size() == y.size() );
	int I = x.size();
	Vec<R> r(I);
	for( int i = 0; i < I; i++ )
		r(i) = x(i) * inv( y(i) );
	return r;
}

/// Matrix-matrix element-wise product
template<class T, class S,
	class R = decltype(typename T::value_type()*typename S::value_type())>
Mat<R> element_prod( const matrix_interface<T> &x, const matrix_interface<S> &y )
{
	assert( x.size1() == y.size1() );
	assert( x.size2() == y.size2() );
	int I = x.size1();
	int J = x.size2();
	Mat<R> r(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			r(i,j) = x(i,j) * y(i,j);
	return r;
};

/// Matrix element-wise inversion
template<class T, class R = typename T::value_type>
Mat<R> element_inv( const matrix_interface<T> &x )
{
	int I = x.size1();
	int J = x.size2();
	Mat<R> r(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			r(i,j) = inv( x(i,j) );
	return r;
}

/// Vector-vector element-wise division
template<class T, class S,
	class R = decltype(typename T::value_type()*typename S::value_type())>
Mat<R> element_div( const matrix_interface<T> &x, const matrix_interface<S> &y )
{
	assert( x.size1() == y.size1() );
	assert( x.size2() == y.size2() );
	int I = x.size1();
	int J = x.size2();
	Mat<R> r(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			r(i,j) = x(i,j) * inv( y(i,j) );
	return r;
};
/// @}

//------------------------------------------------------------------------------

/// @name Transpose and adjoint

/// Transpose of a matrix
//template<class T, class R = typename T::value_type>
//Mat<R> trans( const matrix_interface<T>& A )
//{
//	return ublas::trans(A);
//}
template<class T,
	class R = decltype( ublas::trans( declval<T>() ) )>
matrix_interface<R> trans( const matrix_interface<T>& A )
{
	return ublas::trans(A);
}

/// Adjoint of a vector
template<class T, class R = typename T::value_type>
Vec<R> adjoint( const vector_interface<T>& v )
{
	auto afun = [](const R& x) -> R {return adjoint(x);};
	return entrywise(afun)( v );
}

/// Adjoint of a matrix
template<class T, class R = typename T::value_type>
Mat<R> adjoint( const matrix_interface<T>& A )
{
	auto afun = [](const R& x) -> R {return adjoint(x);};
	return trans( entrywise(afun)( A ) );
}
/// @}

//------------------------------------------------------------------------------
/// @name Inner-product, outer-product and norm

/// Vector inner-product
template<class T, class S,
	class R = decltype(typename T::value_type()*typename S::value_type())>
R inner_prod( const vector_interface<T> &x, const vector_interface<S> &y )
{
	assert(x.size() == y.size());
	int I = x.size();
	R ip  = 0;
	for( int i = 0; i < I; i++ )
		ip += inner_prod( x(i), y(i) );
	return ip;
}

/// Matrix inner-product
template<class T, class S,
	class R = decltype(typename T::value_type()*typename S::value_type())>
R inner_prod( const matrix_interface<T> &x, const matrix_interface<S> &y )
{
	assert(x.size1() == y.size1());
	assert(x.size2() == y.size2());
	int I = x.size1();
	int J = x.size2();
	R r   = 0;
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			r += inner_prod( x(i,j), y(i,j) );
	return r;
}

/// Vector outer-product
//template<class T, class S,
//	class R = decltype(typename T::value_type()*typename S::value_type())>
//Mat<R> outer_prod( const vector_interface<T> &x, const vector_interface<S> &y )
//{
//	return ublas::outer_prod( x, y );
//}
template<class T, class S,
	class R = decltype( ublas::outer_prod( declval<T>(), declval<S>() ) )>
matrix_interface<R> outer_prod( const vector_interface<T> &x, const vector_interface<S> &y )
{
	return ublas::outer_prod( x, y );
}

/// Vector p-norm
template<class T,
	class R = decltype(norm(typename T::value_type())),
	class S = typename T::value_type>
R norm( const vector_interface<T> &x, double p = 2 )
{
	R n;
	auto e_norm = [p]( const S& s ) { return norm( s, p ); };
	if( p == inf )
		if( x.size() == 0 )
			n = 0;
		else
			n = max( entrywise(e_norm)(x) );
	else if( p == 0 )
		n = sum( entrywise(e_norm)(x) );
	else
		n = real( pow( sum( real( pow( entrywise(e_norm)(x), p ) ) ), 1/p ) );
	return n;
}

/// Matrix p-norm
template<class T,
	class R = decltype(norm(typename T::value_type())),
	class S = typename T::value_type>
R norm( const matrix_interface<T> &A, double p = 2 )
{
	R n;
	auto e_norm = [p]( const S& s ) { return norm( s, p ); };
	if( p == inf )
		if( (A.size1() * A.size2()) == 0 )
			n = 0;
		else
			n = max( entrywise(e_norm)(A) );
	else if( p == 0 )
		n = sum( entrywise(e_norm)(A) );
	else
		n = real( pow( sum( real( pow( entrywise(e_norm)(A), p ) ) ), 1/p ) );
	return n;
}
/// @}

/// @}
}
#endif
