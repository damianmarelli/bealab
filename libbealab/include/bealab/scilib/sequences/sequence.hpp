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
class Seq {

	Vec<T> _vec;																/// Vector holding the values
	int _t1;																	/// Initial time of the sequence

public:

	/// @name Constructors
	Seq() : _t1(0) {}
	Seq( int l, int t=0 ) :
		_vec(l), _t1(t) {}
	template<class S>
	Seq( const vector_interface<S>& vec, int t=0 ) :
		_vec(vec), _t1(t) {}
	template<class S>
	Seq( const Seq<S>& seq ) :
		_vec(seq.vec()), _t1(seq.t1()) {}
//	Seq( const initializer_list<T> &l, int t=0 ) :
//		_vec(l), _t1(t) {}
	/// @}

	/// @name Assignment
	template<class I>
	Seq<T>& operator=( const Seq<I> &X )
	{
		_t1  = X.t1();
		_vec = X.vec();
		return *this;
	}
	/// @}

	/// @name Operations
	Seq<T>& operator+=( const Seq<T> &y )
	{
		*this = *this + y;
		return *this;
	}

	Seq<T>& operator+=( const T &a )
	{
		_vec += a;
		return *this;
	}

	Seq<T>& operator-=( const Seq<T> &y )
	{
		*this = *this - y;
		return *this;
	}

	Seq<T>& operator-=( const T &a )
	{
		_vec -= a;
		return *this;
	}

	Seq<T>& operator*=( const Seq<T> &y )
	{
		*this = *this * y;
		return *this;
	}

	Seq<T>& operator*=( const T &a )
	{
		_vec *= a;
		return *this;
	}

//	Seq<T>& operator/=( const Seq<T> &Y )
//	{
//		*this = *this*inv(Y);
//		return *this;
//	}

	Seq<T>& operator/=( const T &a )
	{
		_vec /= a;
		return *this;
	}
	/// @}

	/// @name Get & set
	const Vec<T>& vec() const { return _vec; }
	Vec<T>& vec() { return _vec; }
	int t1() const { return _t1; }
	int t2() const { return _t1 + _vec.size()-1; }

	Seq<T>& t1( int t )
	{
		const_cast<Seq<T>*>(this)->_t1 = t;
		return *const_cast<Seq<T>*>(this);
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
	Seq<T> trunc( int t1n, int t2n ) const
	{
		if( t2n < t1n )
			return Seq<T>();

		// Zero value
		T Z = T();
		if( size() != 0 )
			Z = 0*(*this)(0);

		// Leading zeros
		int N = max( min(_t1, t2n+1) - t1n, 0 );
		Vec<T> v1(N);
		for( int n = 0; n < N; n++ )
			v1(n) = Z;

		// Trailing zeros
		N = max( t2n - max(this->t2(), t1n-1), 0 );
		Vec<T> v3(N);
		for( int n = 0; n < N; n++ )
			v3(n) = Z;

		// Middle part
		Vec<T> v2;
		if( t2n >= _t1 && t1n <= (int)t2() && (int)this->_vec.size() > 0 ) {
			int n1 = min( max( t1n - _t1, 0 ), (int)_vec.size()-1 );
			int n2 = min( max( t2n - _t1, 0 ), (int)_vec.size()-1 );
			v2     = _vec( range( n1, n2+1 ) );
		}

		return Seq<T>( Vec<T>{v1, v2, v3}, t1n );
	}

	Seq<T> trim( double tol )
	{
		if( _vec.size() == 1 && _t1 <= tol )
			return Seq<T>();

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
//	Seq<typename std::result_of<F(T)>::type> apply( F fun ) const
//	{
//		typedef typename std::result_of<F(T)>::type R;
//		int I = this->size();
//		Vec<R> r(I);
//		for( int i = 0; i < I; i++ )
//			r(i) = fun((*this)(i));
//		return {r,t1()};
////		typedef typename std::result_of<F(T)>::type R;
////		const Vec<T>& tmp = dynamic_cast<const Vec<T>&>(*this);
////		return Seq<R>( apply(fun)(tmp), t1() );
//	}

	Seq<T> timereverse() const
	{
		return Seq<T>( flip(_vec), -t2() );
	}

	Seq<T> A() const
	{
		return adjoint(*this);
	}
	/// @}
};

// Types
typedef	Seq<bool>		bseq;
typedef	Seq<int>		iseq;
typedef	Seq<double>		rseq;
typedef	Seq<complex>	cseq;

/// Output
template <class T>
std::ostream &operator<<(std::ostream &os, const Seq<T> &m)
{
	os << "{" << m.vec() << ", " << m.t1() << "}";
	return os;
}

/// Apply a function to all the elements of a sequence
template<class F, class T>
Seq<typename result_of<F(T)>::type>
	entrywise( F fun, const Seq<T>& x )
{
	typedef typename result_of<F(T)>::type R;
	const Vec<T>& x_ = x.vec();
	Vec<R> y_ = entrywise( fun, x_ );
	return { y_, x.t1() };
}

/// Specialization of matvar for Seq<T>
template<class T>
struct matvar<Seq<T>> {

	static
	void* create( const string& varname, const Seq<T>& x )
	{
		Vec<void*> fields(2);
		fields(0) = matvar<double>::create( "t1", x.t1() );
		fields(1) = matvar<Vec<T>>::create( "vec", x.vec() );
		return matfile::create_struct( varname, fields );
	}

	static
	Seq<T> parse( void* pvar )
	{
		void* pt1  = matfile::parse_struct( pvar, "t1" );
		void* pvec = matfile::parse_struct( pvar, "vec" );
		double t1  = matvar<double>::parse( pt1 );
		Vec<T> vec = matvar<Vec<T>>::parse( pvec );
		return Seq<T>( vec, t1 );
	}
};

/// @}
}
#endif
