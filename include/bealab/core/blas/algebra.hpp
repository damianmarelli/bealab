// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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

/// @name Vector operations

/// vector + vector
template< class E1, class E2>
auto operator+( const vector_interface<E1>& x, const vector_interface<E2>& y ) ->
vector_interface<decltype(ublas::operator+(x,y))>
{
	return ublas::operator+( x, y );
}

/// -vector
template<class E>
auto operator-( const vector_interface<E>& x ) ->
vector_interface<decltype(ublas::operator-(x))>
{
	return ublas::operator-( x );
}

/// vector - vector
template<class E1, class E2>
auto operator-( const vector_interface<E1>& x, const vector_interface<E2>& y ) ->
vector_interface<decltype(ublas::operator-(x,y))>
{
	return ublas::operator-( x, y );
}

/// vector * scalar (scalar entries)
template< class E, class T,
	class = typename enable_if< is_scalar<typename E::value_type>::value &&
								is_scalar<T>::value &&
								is_convertible<T,typename E::value_type>::value
							  >::type >
auto operator*( const vector_interface<E>& x, const T& y ) ->
vector_interface<decltype(ublas::operator*(x,y))>
{
	return ublas::operator*( x, y );
}

/// vector * scalar (non-scalar entries)
template< class E, class T,
	class R = decltype( noproxy( declval<typename E::value_type>()*declval<T>() ) ),
	class   = typename enable_if< !is_scalar<typename E::value_type>::value ||
								  !is_scalar<T>::value ||
								  !is_convertible<T,typename E::value_type>::value
								>::type >
Vec<R> operator*( const vector_interface<E>& x, const T& y )
{
	int I = x.size();
	Vec<R> z(I);
	for( int i = 0; i < I; i++ )
		z(i) = x(i) * y;
	return z;
}

/// scalar * vector (scalar entries)
template< class T, class E,
	class = typename enable_if< is_scalar<T>::value &&
								is_scalar<typename E::value_type>::value &&
								is_convertible<T,typename E::value_type>::value
					    	  >::type >
auto operator*( const T& x, const vector_interface<E>& y ) ->
vector_interface<decltype(ublas::operator*(x,y))>
{
	return ublas::operator*( x, y );
}

/// scalar * vector (non-scalar entries)
template< class T, class E,
	class R = decltype( noproxy( declval<T>()*declval<typename E::value_type>() ) ),
	class   = typename enable_if< !is_scalar<T>::value ||
								  !is_scalar<typename E::value_type>::value ||
								  !is_convertible<T,typename E::value_type>::value
								>::type >
Vec<R> operator*( const T& x, const vector_interface<E>& y )
{
	int I = y.size();
	Vec<R> z(I);
	for( int i = 0; i < I; i++ )
		z(i) = x * y(i);
	return z;
}

/// vector / scalar (scalar entries)
template< class E, class T,
	class = typename enable_if< is_scalar<typename E::value_type>::value &&
								is_scalar<T>::value &&
								is_convertible<T,typename E::value_type>::value
							  >::type >
auto operator/( const vector_interface<E>& x, const T& y ) ->
vector_interface<decltype(ublas::operator/(x,y))>
{
	return ublas::operator/( x, y );
}

/// vector / scalar (non-scalar entries)
template< class E, class T,
	class R = decltype( typename E::value_type()/T() ),
	class   = typename enable_if< !is_scalar<T>::value ||
								  !is_scalar<typename E::value_type>::value ||
								  !is_convertible<T,typename E::value_type>::value
								>::type >
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
template< class E1, class E2>
auto operator+( const matrix_interface<E1>& x, const matrix_interface<E2>& y ) ->
matrix_interface<decltype(ublas::operator+(x,y))>
{
	return ublas::operator+(x,y);
}

/// -matrix
template< class E>
auto operator-( const matrix_interface<E>& x ) ->
matrix_interface<decltype(ublas::operator-(x))>
{
	return ublas::operator-(x);
}

/// matrix - matrix
template< class E1, class E2>
auto operator-( const matrix_interface<E1>& x, const matrix_interface<E2>& y ) ->
matrix_interface<decltype(ublas::operator-(x,y))>
{
	return ublas::operator-(x,y);
}

/// matrix * scalar (scalar version)
template< class E, class T,
	class = typename enable_if< is_scalar<typename E::value_type>::value &&
								is_scalar<T>::value &&
								is_convertible<T,typename E::value_type>::value
							  >::type >
auto operator*( const matrix_interface<E>& x, const T& y ) ->
matrix_interface<decltype(ublas::operator*(x,y))>
{
	return ublas::operator*( x, y );
}

/// matrix * scalar (non-scalar version)
template< class E, class T,
	class R = decltype( noproxy( declval<typename E::value_type>()*declval<T>() ) ),
	class   = typename enable_if< !is_scalar<typename E::value_type>::value ||
								  !is_scalar<T>::value ||
								  !is_convertible<T,typename E::value_type>::value
								>::type >
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

/// scalar * matrix (scalar entries)
template< class T, class E,
	class = typename enable_if< is_scalar<T>::value &&
								is_scalar<typename E::value_type>::value &&
								is_convertible<T,typename E::value_type>::value
					    	  >::type >
auto operator*( const T& x, const matrix_interface<E>& y ) ->
matrix_interface<decltype(ublas::operator*(x,y))>
{
	return ublas::operator*( x, y );
}

/// scalar * matrix (non-scalar entries)
template< class T, class E,
	class R = decltype( noproxy( declval<T>()*declval<typename E::value_type>() ) ),
	class   = typename enable_if< !is_scalar<T>::value ||
								  !is_scalar<typename E::value_type>::value ||
								  !is_convertible<T,typename E::value_type>::value
								>::type >
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

/// matrix / scalar (scalar entries)
template< class E, class T,
	class = typename enable_if< is_scalar<typename E::value_type>::value &&
								is_scalar<T>::value &&
								is_convertible<T,typename E::value_type>::value
							  >::type >
auto operator/( const matrix_interface<E>& x, const T& y ) ->
matrix_interface<decltype(ublas::operator/(x,y))>
{
	return ublas::operator/( x, y );
}

/// matrix / scalar (non-scalar entries)
template< class E, class T,
	class R = decltype( typename E::value_type()/T() ),
	class   = typename enable_if< !is_scalar<T>::value ||
								  !is_scalar<typename E::value_type>::value ||
								  !is_convertible<T,typename E::value_type>::value
								>::type >
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

/// matrix * matrix (scalar version)
template< class E1, class E2,
	class = typename enable_if< is_scalar<typename E1::value_type>::value &&
								is_scalar<typename E2::value_type>::value
							  >::type >
auto operator*( const matrix_interface<E1>& x, const matrix_interface<E2>& y ) ->
matrix_interface<decltype(ublas::prod(x,y))>
{
	return ublas::prod( x, y );
}

/// matrix * matrix (non-scalar version)
template< class E1, class E2,
	class R = decltype( noproxy( declval<typename E1::value_type>()*declval<typename E2::value_type>() ) ),
	class   = typename enable_if< !is_scalar<typename E1::value_type>::value ||
								  !is_scalar<typename E2::value_type>::value
								>::type >
Mat<R> operator*( const matrix_interface<E1>& x, const matrix_interface<E2>& y )
{
	assert( x.size2() == y.size1() );
	int I = x.size1();
	int J = y.size2();
	int K = x.size2();
	Mat<R> M(I,J);
	for( int i = 0; i < I; i++ ) {
		for( int j = 0; j < J; j++ ) {
			M(i,j) = R();
			for( int k = 0; k < K; k++ )
				if( k == 0 )
					M(i,j)  = x(i,k) * y(k,j);
				else
					M(i,j) += x(i,k) * y(k,j);
		}
	}
	return M;
}

/// matrix * vector (scalar version)
template< class E1, class E2,
	class = typename enable_if< is_scalar<typename E1::value_type>::value &&
								is_scalar<typename E2::value_type>::value
							  >::type >
auto operator*( const matrix_interface<E1>& x, const vector_interface<E2>& y ) ->
vector_interface<decltype(ublas::prod(x,y))>
{
	return ublas::prod( x, y );
}

/// matrix * vector (non-scalar version)
template< class E1, class E2,
	class R = decltype( noproxy( declval<typename E1::value_type>()*declval<typename E2::value_type>() ) ),
	class   = typename enable_if< !is_scalar<typename E1::value_type>::value ||
								  !is_scalar<typename E2::value_type>::value
								>::type >
Vec<R> operator*( const matrix_interface<E1>& x, const vector_interface<E2>& y )
{
	assert( x.size2() == y.size() );
	int I = x.size1();
	int K = x.size2();
	Vec<R> v(I);
	for( int i = 0; i < I; i++ ) {
		v(i) = R();
		for( int k = 0; k < K; k++ )
			if( k == 0 )
				v(i)  = x(i,k) * y(k);
			else
				v(i) += x(i,k) * y(k);
	}
	return v;
}

/// vector * matrix (scalar version)
template< class E1, class E2,
	class = typename enable_if< is_scalar<typename E1::value_type>::value &&
								is_scalar<typename E2::value_type>::value
							  >::type >
auto operator*( const vector_interface<E1>& x, const matrix_interface<E2>& y ) ->
vector_interface<decltype(ublas::prod(x,y))>
{
	return ublas::prod( x, y );
}

/// vector * matrix (non-scalar version)
template< class E1, class E2,
	class R = decltype( noproxy( declval<typename E1::value_type>()*declval<typename E2::value_type>() ) ),
	class   = typename enable_if< !is_scalar<typename E1::value_type>::value ||
								  !is_scalar<typename E2::value_type>::value
								>::type >
Vec<R> operator*( const vector_interface<E1>& x, const matrix_interface<E2>& y )
{
	assert( x.size() == y.size1() );
	int J = y.size2();
	int K = x.size();
	Vec<R> v(J);
	for( int j = 0; j < J; j++ ) {
		v(j) = R();
		for( int k = 1; k < K; k++ )
			if( k == 0 )
				v(j)  = x(k) * y(k,j);
			else
				v(j) += x(k) * y(k,j);
	}
	return v;
}
/// @}

//------------------------------------------------------------------------------
/// @name Element-wise product/division

/// Vector-vector element-wise product (scalar version)
template< class E1, class E2,
	class = typename enable_if< is_scalar<typename E1::value_type>::value &&
								is_scalar<typename E2::value_type>::value
							  >::type >
auto element_prod( const vector_interface<E1>& x, const vector_interface<E2>& y ) ->
vector_interface<decltype(ublas::element_prod(x,y))>
{
	return ublas::element_prod( x, y );
}

/// Vector-vector element-wise product (non-scalar version)
template< class E1, class E2,
	class R = decltype( noproxy( declval<typename E1::value_type>()*declval<typename E2::value_type>() ) ),
	class   = typename enable_if< !is_scalar<typename E1::value_type>::value ||
								  !is_scalar<typename E2::value_type>::value
								>::type >
Vec<R> element_prod( const vector_interface<E1> &x, const vector_interface<E2> &y )
{
	assert( x.size() == y.size() );
	int I = x.size();
	Vec<R> r(I);
	for( int i = 0; i < I; i++ )
		r(i) = x(i) * y(i);
	return r;
}

/// Vector element-wise inversion
template<class E, class R = typename E::value_type>
Vec<R> element_inv( const vector_interface<E> &x )
{
	int I = x.size();
	Vec<R> r(I);
	for( int i = 0; i < I; i++ )
		r(i) = inv( x(i) );
	return r;
}

/// Vector-vector element-wise division (scalar version)
template< class E1, class E2,
	class = typename enable_if< is_scalar<typename E1::value_type>::value &&
								is_scalar<typename E2::value_type>::value
							  >::type >
auto element_div( const vector_interface<E1>& x, const vector_interface<E2>& y ) ->
vector_interface<decltype(ublas::element_div(x,y))>
{
	return ublas::element_div( x, y );
}

/// Vector-vector element-wise division (non-scalar version)
template< class E1, class E2,
	class R = decltype( noproxy( declval<typename E1::value_type>()/declval<typename E2::value_type>() ) ),
	class   = typename enable_if< !is_scalar<typename E1::value_type>::value ||
								  !is_scalar<typename E2::value_type>::value
								>::type >
Vec<R> element_div( const vector_interface<E1> &x, const vector_interface<E2> &y )
{
	assert( x.size() == y.size() );
	int I = x.size();
	Vec<R> r(I);
	for( int i = 0; i < I; i++ )
		r(i) = x(i) * inv( y(i) );
	return r;
}

/// Matrix-matrix element-wise product (scalar version)
template< class E1, class E2,
	class = typename enable_if< is_scalar<typename E1::value_type>::value &&
								is_scalar<typename E2::value_type>::value
							  >::type >
auto element_prod( const matrix_interface<E1>& x, const matrix_interface<E2>& y ) ->
matrix_interface<decltype(ublas::element_prod(x,y))>
{
	return ublas::element_prod( x, y );
}

/// Matrix-matrix element-wise product (non-scalar version)
template< class E1, class E2,
	class R = decltype( noproxy( declval<typename E1::value_type>()*declval<typename E2::value_type>() ) ),
	class   = typename enable_if< !is_scalar<typename E1::value_type>::value ||
								  !is_scalar<typename E2::value_type>::value
								>::type >
Mat<R> element_prod( const matrix_interface<E1> &x, const matrix_interface<E2> &y )
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
template<class E, class R = typename E::value_type>
Mat<R> element_inv( const matrix_interface<E> &x )
{
	int I = x.size1();
	int J = x.size2();
	Mat<R> r(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			r(i,j) = inv( x(i,j) );
	return r;
}

/// Matrix-matrix element-wise division (scalar version)
template< class E1, class E2,
	class = typename enable_if< is_scalar<typename E1::value_type>::value &&
								is_scalar<typename E2::value_type>::value
							  >::type >
auto element_div( const matrix_interface<E1>& x, const matrix_interface<E2>& y ) ->
matrix_interface<decltype(ublas::element_div(x,y))>
{
	return ublas::element_div( x, y );
}

/// Matrix-matrix  element-wise division (non-scalar version)
template< class E1, class E2,
	class R = decltype( noproxy( declval<typename E1::value_type>()/declval<typename E2::value_type>() ) ),
	class   = typename enable_if< !is_scalar<typename E1::value_type>::value ||
								  !is_scalar<typename E2::value_type>::value
								>::type >
Mat<R> element_div( const matrix_interface<E1> &x, const matrix_interface<E2> &y )
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
template<class T>
auto trans( const matrix_interface<T>& A ) ->
matrix_interface<decltype(ublas::trans(A))>
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
	return entrywise(afun)( trans( A ) );
}
/// @}

//------------------------------------------------------------------------------
/// @name Inner-product, outer-product and norm

/// Vector inner-product (scalar version)
template<class E1, class E2,
	class = typename enable_if< is_scalar<typename E1::value_type>::value &&
								is_scalar<typename E2::value_type>::value
							  >::type >
auto inner_prod( const vector_interface<E1> &x, const vector_interface<E2> &y ) ->
decltype(ublas::inner_prod(x,y))
{
	return ublas::inner_prod(x,y);
}

/// Vector inner-product (non-scalar version)
template<class E1, class E2,
	class R = decltype( noproxy( inner_prod( declval<typename E1::value_type>(), declval<typename E2::value_type>() ) ) ),
	class   = typename enable_if< !is_scalar<typename E1::value_type>::value ||
								  !is_scalar<typename E2::value_type>::value
								>::type >
R inner_prod( const vector_interface<E1> &x, const vector_interface<E2> &y )
{
	assert(x.size() == y.size());
	int I = x.size();
	R ip  = 0;
	for( int i = 0; i < I; i++ )
		ip += inner_prod( x(i), y(i) );
	return ip;
}

/// Matrix inner-product
template<class E1, class E2,
	class R = decltype( noproxy( inner_prod( declval<typename E1::value_type>(),declval<typename E2::value_type>() ) ) ) >
R inner_prod( const matrix_interface<E1> &x, const matrix_interface<E2> &y )
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

/// Vector outer-product (scalar version)
template< class E1, class E2,
	class = typename enable_if< is_scalar<typename E1::value_type>::value &&
								is_scalar<typename E2::value_type>::value
							  >::type >
auto outer_prod( const vector_interface<E1>& x, const vector_interface<E2>& y ) ->
matrix_interface<decltype(ublas::outer_prod(x,y))>
{
	return ublas::outer_prod( x, y );
}

/// Vector outer-product (non-scalar version)
template<class E1, class E2,
	class R = decltype( noproxy( outer_prod( declval<typename E1::value_type>(), declval<typename E2::value_type>() ) ) ),
	class   = typename enable_if< !is_scalar<typename E1::value_type>::value ||
								  !is_scalar<typename E2::value_type>::value
								>::type >
Mat<R> outer_prod( const vector_interface<E1> &x, const vector_interface<E2> &y )
{
	assert(x.size() == y.size());
	int I = x.size();
	Mat<R> op(I,I);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < I; j++ )
			op(i,j) = outer_prod( x(i), y(j) );
	return op;
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
