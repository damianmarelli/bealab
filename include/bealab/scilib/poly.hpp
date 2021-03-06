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
vec<R> poly( const vec<T> &x )
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
vec<R> roots( const vec<T> &x )
{
	int N = x.size() - 1;

	// Companion matrix
	mat<T> A                                = zeros(N,N);
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
R polyval( const vec<T>& p, S x )
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
vec<R> conv( const vec<T>& x, const vec<S>& y )
{
	if( x.size() == 0 || y.size() == 0 )
		return vec<R>();
	int N = x.size() + y.size() - 1;
    vec<R> rv(N);
    for( int n = 0; n < N; n++ )
    	rv(n) = 0*R(x(0)*y(0));
    int Ix = x.size();
    int Iy = y.size();
    for( int ix = 0; ix < Ix; ix++ )
        for( int iy = 0; iy < Iy; iy++ )
            rv(ix+iy) = rv(ix+iy) + R( x(ix) * y(iy) );
    return rv;
}

/// Convolutional matrix
template<class VAL>
mat<VAL> convmat(const vec<VAL> &X, int J)
{
    int L = X.size();
    int I = L + J - 1;
    mat<VAL> M = zeros(I,J);
    vec<VAL> Z = zeros(J-1);
    vec<VAL> B = { Z, flip(X), Z };

    for(int i=0; i<I; i++)
        M.row(i) = B.range( L+J-2-i, L+2*J-3-i );

    return M;
}

/// @}
}
#endif
