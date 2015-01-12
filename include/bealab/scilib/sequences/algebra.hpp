/// @file bealab/scilib/sequences/algebra.hpp
/// Algebraic operations

#ifndef _BEALAB_SEQUENCES_ALGEBRA_
#define	_BEALAB_SEQUENCES_ALGEBRA_

#include <bealab/scilib/poly.hpp>

namespace bealab
{
/// @defgroup sequences_algebra Algebra
/// Algebraic operations
/// @{

/// @name Operations

/// sequence + sequence
template< class T1, class T2,
	class R = decltype( noproxy(T1()+T2()) ) >
Seq<R> operator+( Seq<T1> X, Seq<T2> Y )
{
	if( X.size() == 0 )
		return Y;
	if( Y.size() == 0 )
		return X;

	// Intervals
	int	t1  = min(X.t1(),Y.t1());
	int	ta  = max(X.t1(),Y.t1());
	int	tb  = min(X.t2(),Y.t2());
	int	t2  = max(X.t2(),Y.t2());
	int sx  = X.size();
	int sy  = Y.size();

	// First interval
	Vec<R> v1(0);
	if( t1 == X.t1() && t1 != Y.t1() )
		v1 = X.vec()( range( 0, min(ta-X.t1(),sx) ) );
	else if( t1 != X.t1() && t1 == Y.t1() )
		v1 = Y.vec()( range( 0, min(ta-Y.t1(),sy) ) );

	// Second interval
	Vec<R> v2;
	if( ta <= tb )
		v2 = X.vec()( range( ta-X.t1(), tb-X.t1()+1 ) ) + Y.vec()( range( ta-Y.t1(), tb-Y.t1()+1 ) );
	else {
		int I = ta-tb-1;
		v2.resize(I);
		for( int i = 0; i < I; i++ )
			v2(i) = 0*X.vec()(0);
	}

	// Third interval
	Vec<R> v3(0);
	if( t2 == X.t2() && t2 != Y.t2() )
		v3 = X.vec()( range( max(tb-X.t1()+1,0), X.size() ) );
	else if( t2 != X.t2() && t2 == Y.t2() )
		v3 = Y.vec()( range( max(tb-Y.t1()+1,0), Y.size() ) );

	// Sum
	return Seq<R>( Vec<R>{v1,v2,v3}, t1 );
}

/// -sequence
template< class T >
Seq<T> operator-( const Seq<T> &X )
{
	return Seq<T>( -X.vec(), X.t1() );
}

/// sequence - sequence
template< class T1, class T2,
	class R = decltype( noproxy(T1()-T2()) ) >
Seq<R> operator-(const Seq<T1> &X, const Seq<T2> &Y )
{
	return X + (-Y);
}

/// sequence * sequence
template< class T1, class T2,
	class R = decltype( noproxy(T1()*T2()) ) >
Seq<R> operator*( const Seq<T1> &X, const Seq<T2> &Y )
{
	if( X.size()==0 || Y.size()==0 )
		return Seq<R>();
	Vec<R> vec = conv( X.vec(), Y.vec() );
	int t1     = X.t1() + Y.t1();
	return Seq<R>( vec, t1 );
}

/// sequence * scalar
template< class T1, class T2,
	class R = decltype( noproxy(T1()*T2()) ) >
//	class R = decltype( noproxy(T1()*T2()) ),
//	class = typename enable_if<
//		is_convertible<T2,T1>::value
//	>::type >
Seq<R> operator*( const Seq<T1> &X, const T2 &a )
{
	return Seq<R>( X.vec()*a, X.t1() );
}

/// scalar * sequence
template< class T1, class T2,
	class R = decltype( noproxy(T1()*T2()) ) >
//	class R = decltype( noproxy(T1()*T2()) ),
//	class = typename enable_if<
//		is_convertible<T1,T2>::value
//	>::type >
Seq<R> operator*( const T1 &a, const Seq<T2> &X )
{
	return X*a;
}

/// sequence / scalar
template< class T1, class T2,
	class R = decltype( noproxy(T1()/T2()) ) >
//	class R = decltype( noproxy(T1()/T2()) ),
//	class = typename enable_if<
//		is_convertible<T2,T1>::value
//	>::type >
Seq<R> operator/( const Seq<T1> &X, const T2 &a )
{
    return Seq<R>( X.vec()/a, X.t1() );
}
/// @}

/// @name Element-wise operations
template<class T, class S>
Seq<decltype(T()*S())> element_prod( const Seq<T>& X, const Seq<S>& Y )
{
	typedef decltype(T()*S()) R;

	// If any vector is empty
	if( X.size() == 0 || Y.size() == 0 )
		return Seq<R>();

	// Overlapping interval
	int	ta  = max(X.t1(),Y.t1());
	int	tb  = min(X.t2(),Y.t2());

	// If they do not overlap
	if( ta > tb )
		return Seq<R>();

	// They do overlap
	Vec<R> v = element_prod(
			X.vec()( range( ta-X.t1(), tb-X.t1()+1 ) ),
			Y.vec()( range( ta-Y.t1(), tb-Y.t1()+1 ) ) );
	return Seq<R>( v, ta );
}

template<class T>
Seq<T> element_inv( const Seq<T> &X )
{
    return Seq<T>( element_inv( X.vec() ), X.t1() );
}

template<class T, class S>
Seq<decltype(T()*S())> element_div( const Seq<T>& X, const Seq<S>& Y )
{
	return element_prod( X, element_inv(Y) );
}
/// @}

/// @name Adjoint and Inverse
template<class T>
Seq<T> inv(const Seq<T>& X, int t1, int t2)
{
    int		YL = t2-t1+1;

    // Convolution matrix
    Mat<T>	CX = convmat(X,YL);

    // Impulse function
    int		DL = CX.rows();
    Vec<T>	D = zeros(DL);
    int 	i = -t1-X.t1();
    if(i > 0 && i < DL)
        D(i) = 1;

    // Inverse Seq
    Vec<T>	Y = backslash( CX, D );
    return Seq<T>( Y, t1 );
}

template<class T>
Seq<T> adjoint( Seq<T> X )
{
	int t1 = X.t1();
	int t2 = X.t2();
	for( int t = t1; t <= t2; t++ )
		X(t) = adjoint(X(t));
	return X.timereverse();
}

/// @}

/// @name Inner-product and norm

template<class T>
complex inner_prod(const Seq<T> &X, const Seq<T> &Y )
{
	// If any vector is empty
	if( X.size() == 0 || Y.size() == 0 )
		return 0;

	// Overlapping interval
	int	ta  = max(X.t1(),Y.t1());
	int	tb  = min(X.t2(),Y.t2());

	// If they do not overlap
	if( ta > tb )
		return 0;

	// They do overlap
	complex r = 0;
	for( int t = ta; t <= tb; t++ )
		r += inner_prod( X(t), Y(t) );
	return r;
}

template<class T>
double norm(const Seq<T> &X, double p=2 )
{
	return norm( X.vec(), p );
}
/// @}

/// @}
}
#endif
