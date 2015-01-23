// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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
