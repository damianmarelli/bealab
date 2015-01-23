// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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

/// Convolutional matrix
template<class VAL>
Mat<VAL> convmat(const Vec<VAL> &X, int J)
{
    int L = X.size();
    int I = L + J - 1;
    Mat<VAL> M = zeros(I,J);
    Vec<VAL> Z = zeros(J-1);
    Vec<VAL> B = { Z, flip(X), Z };

    for(int i=0; i<I; i++)
        M.row(i) = B.range( L+J-2-i, L+2*J-3-i );

    return M;
}

/// @}
}
#endif
