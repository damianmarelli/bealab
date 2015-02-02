/// @file bealab/core/blas/comparisons.hpp
/// Operations for entry-wise comparison between vectors and matrices.

#ifndef _BEALAB_BLAS_COMPARISONS_
#define	_BEALAB_BLAS_COMPARISONS_

#include <bealab/core/blas/matrix.hpp>

namespace bealab
{
/// @defgroup blas_comparisons Entry-wise comparisons
/// Operations for entry-wise comparison between vectors and matrices.
/// @{

/// @name Vector-vector comparisons
template<class S, class T>
bvec operator==( const vector_interface<S>& x, const vector_interface<T>& y )
{
	assert( x.size() == y.size() );
	int I = x.size();
	bvec f(I);
	for(int i = 0; i < I; i++)
		f(i) = x(i) == y(i) ? true : false;
	return f;
}

template<class T>
bvec operator==( const vector_interface<T>& x, const typename T::value_type& y )
{
	return x == y*ones(x.size());
}

template<class T>
bvec operator==( const typename T::value_type& y, const vector_interface<T>& x )
{
	return y*ones(x.size()) == x;
}

template<class S, class T>
bvec operator!=( const vector_interface<S>& x, const vector_interface<T>& y )
{
	return ones(x.size())-(x==y);
}

template<class T>
bvec operator!=( const vector_interface<T>& x, const typename T::value_type& y )
{
	return x != y*ones(x.size());
}

template<class T>
bvec operator!=( const typename T::value_type& y, const vector_interface<T>& x )
{
	return y*ones(x.size()) != x;
}

template<class S, class T>
bvec operator<( const vector_interface<S>& x, const vector_interface<T>& y )
{
	assert( x.size() == y.size() );
	int I = x.size();
	bvec f(I);
	for(int i = 0; i < I; i++)
		f(i) = x(i) < y(i) ? true : false;
	return f;
}

template<class T>
bvec operator<( const vector_interface<T>& x, const typename T::value_type& y )
{
	return x < y*ones(x.size());
}

template<class T>
bvec operator<( const typename T::value_type& y, const vector_interface<T>& x )
{
	return y*ones(x.size()) < x;
}

template<class S, class T>
bvec operator>( const vector_interface<S>& x, const vector_interface<T>& y )
{
	return -x < -y;
}

template<class T>
bvec operator>( const vector_interface<T>& x, const typename T::value_type& y )
{
	return x > y*ones(x.size());
}

template<class T>
bvec operator>( const typename T::value_type& y, const vector_interface<T>& x )
{
	return y*ones(x.size()) > x;
}

template<class S, class T>
bvec operator<=( const vector_interface<S>& x, const vector_interface<T>& y )
{
	assert( x.size() == y.size() );
	int I = x.size();
	bvec f(I);
	for(int i = 0; i < I; i++)
		f(i) = x(i) <= y(i) ? true : false;
	return f;
}

template<class T>
bvec operator<=( const vector_interface<T>& x, const typename T::value_type& y )
{
	return x <= y*ones(x.size());
}

template<class T>
bvec operator<=( const typename T::value_type& y, const vector_interface<T>& x )
{
	return y*ones(x.size()) <= x;
}

template<class S, class T>
bvec operator>=( const vector_interface<S>& x, const vector_interface<T>& y )
{
	return -x <= -y;
}

template<class T>
bvec operator>=( const vector_interface<T>& x, const typename T::value_type& y )
{
	return x >= y*ones(x.size());
}

template<class T>
bvec operator>=( const typename T::value_type& y, const vector_interface<T>& x )
{
	return y*ones(x.size()) >= x;
}
/// @}

/// @}
}
#endif
