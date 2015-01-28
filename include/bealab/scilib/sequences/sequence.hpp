// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/scilib/sequences/sequence.hpp
/// Defines a templated class of sequences with compact support.

#ifndef _BEALAB_SEQUENCES_SEQUENCE_
#define	_BEALAB_SEQUENCES_SEQUENCE_

#include <bealab/core/blas.hpp>
#include <bealab/core/matfile.hpp>

namespace bealab
{
/// @defgroup sequences_class Sequences
/// Defines a templated class of sequences with compact support.
/// @{

/// Sequence class
template<class T>
class sequence {

	vec<T> _vec;																/// Vector holding the values
	int _t1;																	/// Initial time of the sequence

public:

	/// @name Constructors
	sequence() : _t1(0) {}
	sequence( int l, int t=0 ) :
		_vec(l), _t1(t) {}
	template<class S>
	sequence( const vector_interface<S>& vec, int t=0 ) :
		_vec(vec), _t1(t) {}
	template<class S>
	sequence( const sequence<S>& seq ) :
		_vec(seq.buffer()), _t1(seq.t1()) {}
//	sequence( const initializer_list<T> &l, int t=0 ) :
//		_vec(l), _t1(t) {}
	/// @}

	/// @name Assignment
	template<class I>
	sequence<T>& operator=( const sequence<I> &X )
	{
		_t1  = X.t1();
		_vec = X.buffer();
		return *this;
	}
	/// @}

	/// @name Operations
	sequence<T>& operator+=( const sequence<T> &y )
	{
		*this = *this + y;
		return *this;
	}

	sequence<T>& operator+=( const T &a )
	{
		_vec += a;
		return *this;
	}

	sequence<T>& operator-=( const sequence<T> &y )
	{
		*this = *this - y;
		return *this;
	}

	sequence<T>& operator-=( const T &a )
	{
		_vec -= a;
		return *this;
	}

	sequence<T>& operator*=( const sequence<T> &y )
	{
		*this = *this * y;
		return *this;
	}

	sequence<T>& operator*=( const T &a )
	{
		_vec *= a;
		return *this;
	}

//	sequence<T>& operator/=( const sequence<T> &Y )
//	{
//		*this = *this*inv(Y);
//		return *this;
//	}

	sequence<T>& operator/=( const T &a )
	{
		_vec /= a;
		return *this;
	}
	/// @}

	/// @name Get & set
	const vec<T>& buffer() const { return _vec; }
	vec<T>& buffer() { return _vec; }
	int t1() const { return _t1; }
	int t2() const { return _t1 + _vec.size()-1; }

	sequence<T>& t1( int t )
	{
		const_cast<sequence<T>*>(this)->_t1 = t;
		return *const_cast<sequence<T>*>(this);
	}

	int  size() const { return _vec.size(); }

	T& operator()( int t )
	{
		int _t2 = t2();
		if( t < _t1 || t > _t2 )
			*this = trunc( min(_t1,t), max(_t2,t) );
		return _vec( t - _t1 );
	}

	T operator()( int t ) const
	{
		int _t2 = t2();
		if( t < _t1 || t > _t2 )
			return T();
		else
			return _vec( t - _t1 );
	}
	/// @}

	/// @name Other operations
	sequence<T> trunc( int t1n, int t2n ) const
	{
		if( t2n < t1n )
			return sequence<T>();

		// Zero value
		T Z = T();
		if( size() != 0 )
			Z = 0*(*this)(0);

		// Leading zeros
		int N = max( min(_t1, t2n+1) - t1n, 0 );
		vec<T> v1(N);
		for( int n = 0; n < N; n++ )
			v1(n) = Z;

		// Trailing zeros
		N = max( t2n - max(this->t2(), t1n-1), 0 );
		vec<T> v3(N);
		for( int n = 0; n < N; n++ )
			v3(n) = Z;

		// Middle part
		vec<T> v2;
		if( t2n >= _t1 && t1n <= (int)t2() && (int)this->_vec.size() > 0 ) {
			int n1 = min( max( t1n - _t1, 0 ), (int)_vec.size()-1 );
			int n2 = min( max( t2n - _t1, 0 ), (int)_vec.size()-1 );
			v2     = _vec( range( n1, n2+1 ) );
		}

		return sequence<T>( vec<T>{v1, v2, v3}, t1n );
	}

	sequence<T> trim( double tol )
	{
		if( _vec.size() == 1 && _t1 <= tol )
			return sequence<T>();

		int t1 = this->t1();
		int t2 = this->t2();
		int ta = t1;
		int tb = t2;
		for( int t = t1; t <= t2; t++ ) {
			ta = t;
			if( norm((*this)(t)) >= tol )
				break;
		}
		for( int t = t2; t >= t1; t-- ) {
			tb = t;
			if( norm((*this)(t)) >= tol )
				break;
		}
		return trunc( ta, tb );
	}

//	template<class F>
//	sequence<typename std::result_of<F(T)>::type> apply( F fun ) const
//	{
//		typedef typename std::result_of<F(T)>::type R;
//		int I = this->size();
//		vec<R> r(I);
//		for( int i = 0; i < I; i++ )
//			r(i) = fun((*this)(i));
//		return {r,t1()};
////		typedef typename std::result_of<F(T)>::type R;
////		const vec<T>& tmp = dynamic_cast<const vec<T>&>(*this);
////		return sequence<R>( apply(fun)(tmp), t1() );
//	}

	sequence<T> timereverse() const
	{
		return sequence<T>( flip(_vec), -t2() );
	}

	sequence<T> A() const
	{
		return adjoint(*this);
	}
	/// @}
};

// Types
typedef	sequence<bool>		bseq;
typedef	sequence<int>		iseq;
typedef	sequence<double>	rseq;
typedef	sequence<complex>	cseq;

/// Output
template <class T>
std::ostream &operator<<(std::ostream &os, const sequence<T> &m)
{
	os << "{" << m.buffer() << ", " << m.t1() << "}";
	return os;
}

/// Apply a function to all the elements of a sequence
template<class F, class T>
sequence<typename result_of<F(T)>::type>
	entrywise( F fun, const sequence<T>& x )
{
	typedef typename result_of<F(T)>::type R;
	const vec<T>& x_ = x.buffer();
	vec<R> y_ = entrywise( fun, x_ );
	return { y_, x.t1() };
}

/// Specialization of matvar for sequence<T>
template<class T>
struct matvar<sequence<T>> {

	static
	void* create( const string& varname, const sequence<T>& x )
	{
		vec<void*> fields(2);
		fields(0) = matvar<double>::create( "t1", x.t1() );
		fields(1) = matvar<vec<T>>::create( "vec", x.buffer() );
		return matfile::create_struct( varname, fields );
	}

	static
	sequence<T> parse( void* pvar )
	{
		void* pt1 = matfile::parse_struct( pvar, "t1" );
		void* pv  = matfile::parse_struct( pvar, "vec" );
		double t1 = matvar<double>::parse( pt1 );
		vec<T> v  = matvar<vec<T>>::parse( pv );
		return sequence<T>( v, t1 );
	}
};

/// @}
}
#endif
