// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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
sequence<R> operator+( sequence<T1> X, sequence<T2> Y )
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
	vec<R> v1(0);
	if( t1 == X.t1() && t1 != Y.t1() )
		v1 = X.buffer()( range( 0, min(ta-X.t1(),sx) ) );
	else if( t1 != X.t1() && t1 == Y.t1() )
		v1 = Y.buffer()( range( 0, min(ta-Y.t1(),sy) ) );

	// Second interval
	vec<R> v2;
	if( ta <= tb )
		v2 = X.buffer()( range( ta-X.t1(), tb-X.t1()+1 ) ) + Y.buffer()( range( ta-Y.t1(), tb-Y.t1()+1 ) );
	else {
		int I = ta-tb-1;
		v2.resize(I);
		for( int i = 0; i < I; i++ )
			v2(i) = 0*X.buffer()(0);
	}

	// Third interval
	vec<R> v3(0);
	if( t2 == X.t2() && t2 != Y.t2() )
		v3 = X.buffer()( range( max(tb-X.t1()+1,0), X.size() ) );
	else if( t2 != X.t2() && t2 == Y.t2() )
		v3 = Y.buffer()( range( max(tb-Y.t1()+1,0), Y.size() ) );

	// Sum
	return sequence<R>( vec<R>{v1,v2,v3}, t1 );
}

/// -sequence
template< class T >
sequence<T> operator-( const sequence<T> &X )
{
	return sequence<T>( -X.buffer(), X.t1() );
}

/// sequence - sequence
template< class T1, class T2,
	class R = decltype( noproxy(T1()-T2()) ) >
sequence<R> operator-(const sequence<T1> &X, const sequence<T2> &Y )
{
	return X + (-Y);
}

/// sequence * sequence
template< class T1, class T2,
	class R = decltype( noproxy(T1()*T2()) ) >
sequence<R> operator*( const sequence<T1> &X, const sequence<T2> &Y )
{
	if( X.size()==0 || Y.size()==0 )
		return sequence<R>();
	vec<R> vec = conv( X.buffer(), Y.buffer() );
	int t1     = X.t1() + Y.t1();
	return sequence<R>( vec, t1 );
}

/// sequence * scalar
template< class T1, class T2,
	class R = decltype( noproxy(T1()*T2()) ) >
//	class R = decltype( noproxy(T1()*T2()) ),
//	class = typename enable_if<
//		is_convertible<T2,T1>::value
//	>::type >
sequence<R> operator*( const sequence<T1> &X, const T2 &a )
{
	return sequence<R>( X.buffer()*a, X.t1() );
}

/// scalar * sequence
template< class T1, class T2,
	class R = decltype( noproxy(T1()*T2()) ) >
//	class R = decltype( noproxy(T1()*T2()) ),
//	class = typename enable_if<
//		is_convertible<T1,T2>::value
//	>::type >
sequence<R> operator*( const T1 &a, const sequence<T2> &X )
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
sequence<R> operator/( const sequence<T1> &X, const T2 &a )
{
    return sequence<R>( X.buffer()/a, X.t1() );
}
/// @}

/// @name Element-wise operations
template<class T, class S>
sequence<decltype(T()*S())> element_prod( const sequence<T>& X, const sequence<S>& Y )
{
	typedef decltype(T()*S()) R;

	// If any vector is empty
	if( X.size() == 0 || Y.size() == 0 )
		return sequence<R>();

	// Overlapping interval
	int	ta  = max(X.t1(),Y.t1());
	int	tb  = min(X.t2(),Y.t2());

	// If they do not overlap
	if( ta > tb )
		return sequence<R>();

	// They do overlap
	vec<R> v = element_prod(
			X.buffer()( range( ta-X.t1(), tb-X.t1()+1 ) ),
			Y.buffer()( range( ta-Y.t1(), tb-Y.t1()+1 ) ) );
	return sequence<R>( v, ta );
}

template<class T>
sequence<T> element_inv( const sequence<T> &X )
{
    return sequence<T>( element_inv( X.buffer() ), X.t1() );
}

template<class T, class S>
sequence<decltype(T()*S())> element_div( const sequence<T>& X, const sequence<S>& Y )
{
	return element_prod( X, element_inv(Y) );
}
/// @}

/// @name Adjoint and Inverse
template<class T>
sequence<T> inv(const sequence<T>& X, int t1, int t2)
{
    int		YL = t2-t1+1;

    // Convolution matrix
    mat<T>	CX = convmat(X,YL);

    // Impulse function
    int		DL = CX.rows();
    vec<T>	D = zeros(DL);
    int 	i = -t1-X.t1();
    if(i > 0 && i < DL)
        D(i) = 1;

    // Inverse sequence
    vec<T>	Y = backslash( CX, D );
    return sequence<T>( Y, t1 );
}

template<class T>
sequence<T> adjoint( sequence<T> X )
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
complex inner_prod(const sequence<T> &X, const sequence<T> &Y )
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
double norm(const sequence<T> &X, double p=2 )
{
	return norm( X.buffer(), p );
}
/// @}

/// @}
}
#endif
