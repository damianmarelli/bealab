/// @file bealab/scilib/poly.hpp
/// Polynomimals

#ifndef _BEALAB_POLY_
#define	_BEALAB_POLY_

#include <bealab/core/blas/vector.hpp>
#include <bealab/core/blas/matrix.hpp>

namespace bealab
{
/// @defgroup poly Polynomimals
/// Polynomimal operations. A polynomial 'p(x)' is represented by a vector 'a'
/// such that p(x) = a(0) x^n + a(1) x^(n-1) + ... + a(n).
/// @{

/// Obtain a polynomial given its roots. The call a = poly(r) returns the
/// coefficients of the polynomial p(x) = a(0) x^n + a(1) x^(n-1) + ... + a(n),
/// whose roots are in 'r'.
template<class R=complex,class T>
Vec<R> poly( const Vec<T> &x )
{
	int N  = x.size();
	cvec c = zeros(N+1);
	c(0)   = 1;
	for( int j = 0; j < N; j++ )
		c( range(1,j+2) ) = c( range(1,j+2) ) - x(j) * c( range(0,j+1) );
	return c;
}

/// Obtain the roots of a a polynomial. The roots of the polynomial
/// p(x) = a(0) x^n + a(1) x^(n-1) + ... + a(n) are obtained by r = poly(a)
template<class R=complex,class T>
Vec<R> roots( const Vec<T> &x )
{
	int N = x.size() - 1;

	// Companion matrix
	Mat<T> A                                = zeros(N,N);
	A.range( 1, A.size1()-1, 0, A.size2()-2 ) = eye(N-1);
	A.col( A.size2()-1 )                    = - flip( x.range(1,x.size()-1) ) / x(0);

	// Roots
	cmat D = get<0>( eig(A) );
	return diag( D );
}

/// Evaluate a polynomial at a given value. The call y = poly( a, x )
/// returns y = a(0) x^n + a(1) x^(n-1) + ... + a(n).
template<class T, class S,
	class R = decltype(T()*S()) >
R polyval( const Vec<T>& p, S x )
{
	int N = p.size();
	R rv  = 0;
	for( int n = 0; n < N; n++ )
		rv += p(n) * pow(x,N-n-1);
	return rv;
}

/// Polynomial multiplication
template<class T, class S,
	class R = decltype( noproxy(T()*S()) ) >
Vec<R> conv( const Vec<T>& x, const Vec<S>& y )
{
	if( x.size() == 0 || y.size() == 0 )
		return Vec<R>();
	int N = x.size() + y.size() - 1;
    Vec<R> rv(N);
    for( int n = 0; n < N; n++ )
    	rv(n) = 0*R(x(0)*y(0));
    int Ix = x.size();
    int Iy = y.size();
    for( int ix = 0; ix < Ix; ix++ )
        for( int iy = 0; iy < Iy; iy++ )
            rv(ix+iy) = rv(ix+iy) + R( x(ix) * y(iy) );
    return rv;
}

//// XXX esta delcaracion se necesita aca abajo
//template<class T>
//cvec ifft( const ublas::vector_expression<T>& x );
//
//template<class T, class S>
//auto convf( const ublas::vector_expression<T>& xe, const ublas::vector_expression<S>& ye ) ->
//Vec<decltype(xe()(0)*ye()(0))>
//{
//	// XXX eliminar estos casts cuando la ifft accepte un vector_expression
//	Vec<typename T::value_type> x = xe();
//	Vec<typename T::value_type> y = ye();
//	if( x.size() == 0 || y.size() == 0 )
//		error("conv() - vector with zero size");
//	typedef decltype(xe()(0)*ye()(0)) R;
//	int N = x.size() + y.size() - 1;
//
//	const cvec& xf = dtft( x, N );
//	const cvec& yf = dtft( y, N );
//	const cvec& zf = element_prod( xf, yf );
//	return ifft( zf );
//}

//template<class T, class S>
//auto conv( const ublas::vector_expression<T>& xe, const ublas::vector_expression<S>& ye ) ->
//Vec<decltype(xe()(0)*ye()(0))>
//{
//	const int conv_switch
//}

/// @}
}
#endif
