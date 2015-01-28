// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/core/blas/vector.hpp
/// Vector modeling

#ifndef _BEALAB_BLAS_VECTOR_
#define	_BEALAB_BLAS_VECTOR_

#include <bealab/core/blas/prelim.hpp>

namespace bealab
{
/// @defgroup blas_vector Vectors
/// Vector modeling
/// @{

/// Vector class with extended construction
template<class value_type>
class dense_vector : public ublas::vector<value_type> {
public:

	using ublas::vector<value_type>::vector;
	using ublas::vector<value_type>::operator=;

	/// Default constructor
	dense_vector() = default;

	/// Construct a vector with the data between two given iterators
	template<class It>
	dense_vector( It first, It last )
	{
		this->resize( distance( first, last ) );
		int n = 0;
		for( It it = first; it < last; it++, n++ )
			(*this)(n) = *it;
	}

	/// Construct a vector with the data given in an initializer_list
	dense_vector( const initializer_list<value_type>& list )
	{
		this->resize( list.size() );
		int i = 0;
		for( auto it = list.begin(); it != list.end(); it++ )
			(*this)(i++) = *it;
	}

	/// Construct a vector by concatenating the sub-vectors given in an
	/// initializer_list
	dense_vector( const initializer_list<dense_vector<value_type>> &list )
	{
		// Init
		auto p = list.begin();
		int I  = list.size();

		// Compute total size
		int L = 0;
	    for(int i = 0; i < I; i++)
	        L += p[i].size();

	    // Fill the vector
	    this->resize(L);
	    int cursor = 0;
	    for(int i = 0; i < I; i++) {
	    	int l = p[i].size();
	    	ublas::vector_range<ublas::vector<value_type>>( *this,
	    			range(cursor,cursor+l) ) = p[i];
	    	cursor += l;
	    }
	}

	/// Copy into the vector the data given in an initializer_list
	dense_vector<value_type>& operator=( const initializer_list<value_type>& list )
	{
		*this = dense_vector<value_type>(list);
		return *this;
	}

	/// Copy into the vector the data obtained by concatenating the sub-vectors
	/// given in an initializer_list
	dense_vector<value_type>& operator=( const initializer_list<dense_vector<value_type>> &list )
	{
		*this = dense_vector<value_type>(list);
		return *this;
	}
};

/// Interface for any vector_expression
template<class base>
class vector_interface : public base {
public:

	using typename base::value_type;
	using typename base::const_reference;
	using typename base::reference;
	using base::base;
	using base::operator=;

	/// Default constructor (because it is not inherited)
	vector_interface() = default;

	/// Copy constructor using the base
	vector_interface( const base& b ) : base(b) {}

	/// Access the base (constant)
	const base& operator()() const
	{
		return *this;
	}

	/// Access the base (non-constant)
	base& operator()()
	{
		return *this;
	}

	/// Access an entry (constant)
	const_reference operator()( int n ) const
	{
		return base::operator()( n );
	}

	/// Access an entry (non-constant)
	reference operator()( int n )
	{
		return base::operator()( n );
	}

	/// Access the entries in a range
	vector_interface<ublas::vector_range<base>>
	operator()( const ublas::range& r ) const
	{
		auto& self = const_cast<vector_interface<base>&>(*this);
		return { self, r };
	}

	/// Access the entries in a slice
	vector_interface<ublas::vector_slice<base>>
	operator()( const ublas::slice& s ) const
	{
		auto& self = const_cast<vector_interface<base>&>(*this);
		return { self, s };
	}

	/// Indirect entry access
	vector_interface<ublas::vector_indirect<base>>
	operator()( const indirect& idxs ) const
	{
		auto& self = const_cast<vector_interface<base>&>(*this);
		return { self, idxs };
	}
};

/// Dense vector template
template<class value_type>
using vec = vector_interface<dense_vector<value_type>>;

/// Sparse vector template
template<class value_type>
using sparse_vector = vector_interface<ublas::compressed_vector<value_type>>;

// Shorthand expressions for dense vectors
typedef	vec<bool> bvec;															///< Boolean vector
typedef	vec<int> ivec;															///< Integer vector
typedef vec<double> rvec;														///< Real vector
typedef	vec<complex> cvec;														///< Complex vector

/// Convert an expression template into a temporary
template<class E>
vector_interface<dense_vector<typename E::value_type>>
noproxy( const vector_interface<E>& x )
{
	return x;
}

/// Displays a vector in the console
template<class A, class B, class T>
std::basic_ostream<A,B> &operator<<( std::basic_ostream<A,B> &os,
		                             const vector_interface<T> &v )
{
    int I = v.size();
    os << '[';
    for( int i = 0; i < I-1; i++ )
        os << v(i) << ", ";
    if( I > 0 )
    	os << v(I-1);
    os << ']';
    return os;
}

/// Apply a function to all the elements of a vector
template<class F, class T>
vec<typename result_of<F(typename T::value_type)>::type>
	entrywise( F fun, const vector_interface<T>& x )
{
	typedef typename result_of<F(typename T::value_type)>::type R;
	int N = x.size();
	vec<R> y(N);
	for( int n = 0; n < N; n++ )
		y(n) = fun( x(n) );
	return y;
}

/// @}
}
#endif
