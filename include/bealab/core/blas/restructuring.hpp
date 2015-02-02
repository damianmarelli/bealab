/// @file bealab/core/blas/restructuring.hpp
/// Operations for restructuring vectors and matrices.

#ifndef _BEALAB_BLAS_RESTRUCTURING_
#define	_BEALAB_BLAS_RESTRUCTURING_

#include <bealab/core/blas/matrix.hpp>

namespace bealab
{
/// @defgroup blas_restructuring Restructuring
/// Operations for restructuring vectors and matrices.
/// @{

/// Make a vector with the diagonal of a matrix
template<class T>
vector_interface<ublas::matrix_vector_range<T>> diag( const matrix_interface<T>& A )
{
	int N       = min( A.size1(), A.size2() );
	const T& A_ = dynamic_cast<const T&>( A );
	return ublas::matrix_vector_range<T>( const_cast<T&>( A_ ), range(0,N), range(0,N) );
}

/// Make a diagonal matrix with the elements of a vector
template<class T>
banded_matrix<typename T::value_type>
diag( const vector_interface<T>& x )
{
	typedef typename T::value_type R;
	int N = x.size();
	banded_matrix<R> rv(N,N,0,0);
	for( int n = 0; n < N; n++ )
		rv(n,n) = x(n);
	return rv;
}

/// Flip the elements of a vector
template<class T>
vector_interface<ublas::vector_indirect<T>>
flip( const vector_interface<T> &x )
{
	ivec idxs   = vslice( (int)x.size()-1, -1, (int)x.size() );
	const T& x_ = dynamic_cast<const T&>( x );
	return { const_cast<T&>(x_), indirect(idxs) };
}

/// Circular-shift of the elements of a vector.
/// The shift is done on the right by 'n' elements.
template<class T>
vector_interface<ublas::vector_indirect<T>>
circshift( const vector_interface<T> &x, int n )
{
	int N = x.size();
	int s = mod( n, N );
	ivec idxs = { vrange( N-s, N ), vrange( 0, N-s ) };
	const T& x_ = dynamic_cast<const T&>( x );
	return { const_cast<T&>(x_), indirect(idxs) };
}

/// Circular-shift of the elements of a matrix.
/// The shift is done down by 'm' elements and on the right by 'n' elements.
template<class T>
matrix_interface<ublas::matrix_indirect<T>>
circshift( const matrix_interface<T> &A, int m, int n )
{
	int M      = A.size1();
	int N      = A.size2();
	int sm     = mod( m, M );
	int sn     = mod( n, N );
	ivec idxsm = { vrange( M-sm, M ), vrange( 0, M-sm ) };
	ivec idxsn = { vrange( N-sn, N ), vrange( 0, N-sn ) };
	return { const_cast<T&>(A), indirect(idxsm), indirect(idxsn) };
}

/// @}
}
#endif
