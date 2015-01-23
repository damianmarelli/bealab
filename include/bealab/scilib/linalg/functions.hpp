// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/scilib/linalg/functions.hpp
/// Functions for linear algebra

#ifndef _BEALAB_LINALG_FUNCTIONS_
#define	_BEALAB_LINALG_FUNCTIONS_

#include <bealab/core/blas.hpp>
#include <bealab/scilib/linalg/lapack.hpp>

namespace bealab
{
/// @defgroup linalg_functions Functions
/// Functions for linear algebra
/// @{

/// @name Equation solving

/// Solve LLS problem.
/// Call: X = lls( A, Y ). The returned matrix X is the minimum norm solution
/// of || Y - A * X ||
template<class E, class F,
	class R = decltype( noproxy( declval<typename E::value_type>()*declval<typename F::value_type>() ) ) >
Mat<R> lls( const matrix_interface<E>& A, const matrix_interface<F>& Y )
{
	// Assert dimensions
	assert( A.size1() == Y.size1() );

	// Deal with pathological cases
	if( Y.size1() == 0 || Y.size2() == 0 || A.size2() == 0 )
		return zeros( A.size2(), Y.size2() );

	// Convert to column major storage
	ublas::matrix<R,ublas::column_major> A_ = A;
	ublas::matrix<R,ublas::column_major> Y_ = Y;

	// Call parameters
	int M    = A_.size1();
	int N    = A_.size2();
	int Nrhs = Y_.size2();
	if( M == 0 )
		A_ = zeros(1,N);
	int lda  = max( 1, M );
	int ldb  = max( 1, max(M,N) );
	if( ldb > M )
		Y_.resize( ldb, Nrhs );
	int minmn = min(M,N);
	rvec s(minmn);
	double rcond = -1;
	int rank;
	int smlsiz = lapack::ilaenv<R>( 9, "GELSD", "", 0, 0, 0, 0 );
	int nlvl = max( 0, int( log2( min( M,N )/(smlsiz+1) ) ) + 1 );
	int liwork = 3 * minmn * nlvl + 11 * minmn;
	int iwork[max(1,liwork)];

	// Ask for the optimal value of lwork
	int lwork = -1;
	R qwork;
	if( lapack::gelsd( M, N, Nrhs, &A_(0,0), lda, &Y_(0,0), ldb, &s(0), rcond,
			rank, &qwork, lwork, iwork ) )
		error("lls() - Error calling lapack::gelsd()");

	// Solve the system
	lwork = real(qwork);
	R work[max(1,lwork)];
	if( lapack::gelsd( M, N, Nrhs, &A_(0,0), lda, &Y_(0,0), ldb, &s(0), rcond,
			rank, work, lwork, iwork ) )
		error("lls() - Error calling lapack::gelsd()");

	// Parse the solution
	return ublas::matrix_range<decltype(Y_)>( Y_, range(0,N), range(0,Nrhs) );
}

/// Solve LLS problem.
/// Call: x = lls( A, y ). The returned vector x is the minimum norm solution
/// of || y - A * x ||
template<class E, class F,
	class R = decltype( noproxy( declval<typename E::value_type>()*declval<typename F::value_type>() ) ) >
Vec<R> lls( const matrix_interface<E>& A, const vector_interface<F>& y )
{
	typedef typename F::value_type U;
	Mat<U> Y( y.size(), 1 );
	Y.column(0) = y;
	return lls( A, Y ).column(0);
}

/// Solve linear system of equations.
/// Call: X = linsolve( A, Y ). Returns a matrix X such that Y = A * X.
template<class E, class F,
	class R = decltype( noproxy( declval<typename E::value_type>()*declval<typename F::value_type>() ) ) >
Mat<R> linsolve( const matrix_interface<E>& A, const matrix_interface<F>& Y )
{
	// Convert to column major storage
	ublas::matrix<R,ublas::column_major> A_ = A;
	ublas::matrix<R,ublas::column_major> Y_ = Y;

	// Assert dimensions
	assert( A_.size1() == A_.size2() );
	assert( A_.size2() == Y_.size1() );

	// Deal with pathological cases
	if( A_.size1() == 0 || Y_.size2() == 0 )
		return zeros( A_.size2(), Y_.size2() );

	// Call parameters
	int N    = A_.size1();
	int Nrhs = Y_.size2();
	int lda  = max( 1, N );
	int ipiv[N];
	int ldb  = max( 1, N );
	if( lapack::gesv( N, Nrhs, &A_(0,0), lda, ipiv, &Y_(0,0), ldb ) )
		error("linsolve() - Error calling lapack::gesv()");

	return Y_;
}

/// Solve linear system of equations.
/// Call: x = lsolve( A, b ). Returns a vector x such that y = A * x.
template<class E, class F,
	class R = decltype( noproxy( declval<typename E::value_type>()*declval<typename F::value_type>() ) ) >
Vec<R> linsolve( const matrix_interface<E>& A, const vector_interface<F>& x )
{
	typedef typename F::value_type U;
	Mat<U> X( x.size(), 1 );
	X.column(0) = x;
	return linsolve( A, X ).column(0);
}
/// @}

/// @name Decompositions

/// SVD decomposition.
/// Call: USV = svd(A). Returns a tuple USV with U, S and V such that A = U * S * V.A()
template<class E, class T = typename E::value_type>
tuple< Mat<T>, rmat, Mat<T> > svd( const matrix_interface<E>& A )
{
	// Convert to column major storage
	ublas::matrix<T,ublas::column_major> A_ = A;

	// Call parameters
	int M = A_.size1();
	int N = A_.size2();
	ublas::matrix<T,ublas::column_major> U_(M,M), VT_(N,N);
	rvec s(min(N,M));
	char jobz = 'A';
	int lda   = max( 1, M );
	int ldu   = M;
	int lvdt  = N;
	int iwork[8*min(M,N)];

	// Ask for the optimal value of lwork
	int lwork = -1;
	T qwork;
	if( lapack::gesdd( jobz, M, N, &A_(0,0), lda, &s(0), &U_(0,0), ldu,
			 &VT_(0,0), lvdt, &qwork, lwork, iwork ) )
		error("svd() - Error calling lapack::gesdd()");

	// Compute the SVD
	lwork = real(qwork);
	T work[lwork];
	if( lapack::gesdd( jobz, M, N, &A_(0,0), lda, &s(0), &U_(0,0), ldu,
			 &VT_(0,0), lvdt, work, lwork, iwork ) )
		error("svd() - Error calling lapack::gesdd()");

	// Parse results
	Mat<T> U = U_;
	Mat<T> V = adjoint(Mat<T>(VT_));
	rmat S   = diag( s );
	if( M > N )
		S = rmat( {{S}, {zeros(M-N,N)}});
	else if( M < N )
		S = rmat( {{S, zeros(M,N-M)}} );
	return tuple< Mat<T>, rmat, Mat<T> >{ U, S, V };
}

/// Eigenvalue decomposition.
/// Call: DT = eig(A). Returns a tuple DT with D and T such that A = T * D * inv(T)
template<class E>
tuple< cmat, cmat > eig( const matrix_interface<E>& A )
{
	// Convert to column major storage
	typedef complex T;		// XXX Patch because the LAPACK function for real matrices doesn't work
	ublas::matrix<T,ublas::column_major> A_ = A;

	// Assert dimensions
	assert( A_.size1() == A_.size2() );

	// Deal with the zero-size case
	if( A_.size1() == 0 )
		return tuple<cmat,cmat>{ cmat(), cmat() };

	// Call parameters
	char jobl = 'N';
	char jobr = 'V';
	int N     = A_.size1();
	int lda   = max(1,N);
	cvec W(N);
	int ldvl = 1;
	cmat VL(ldvl,N);
	int ldvr = N;
	cmat VR(ldvr,N);

	// Ask for the optimal value of lwork
	int lwork = -1;
	T qwork;
	if( lapack::geev( jobl, jobr, N, &A_(0,0), lda, &W(0), &VL(0,0), ldvl,
			&VR(0,0), ldvr, &qwork, lwork ) )
		error("eig() - Error calling lapack::geev()");

	// Lapack call
	lwork = real(qwork);
	T work[max(1,lwork)];
	if( lapack::geev( jobl, jobr, N, &A_(0,0), lda, &W(0), &VL(0,0), ldvl,
			&VR(0,0), ldvr, work, lwork ) )
		error("eig() - Error calling lapack::geev()");

	// Parse
	return tuple<cmat,cmat>{ diag(W), trans(VR) };
}

/// Cholesky factorization.
/// Call X = chol( A ). Returns an upper-triangular matrix X such that A = X.A() * X.
/// A has to be positive. It assumes that A is Hermitian, so it only considers
/// its diagonal and upper triangle.
template<class E, class T = typename E::value_type>
Mat<T> cholesky( const matrix_interface<E>& A )
{
	typedef ublas::matrix<T,ublas::column_major> CM;
	CM A_ = A;

	// Assert dimensions and self-adjointness
	assert( A.size1() == A.size2() );
//#ifndef NDEBUG
//	assert( norm( A-A.A() ) == 0 );
//#endif

	// Call parameters
	char uplo = 'U';
	int N     = A_.size1();
	int lda   = max( 1, N );
	int info  = lapack::potrf( uplo, N, &A_(0,0), lda );
	if( info ) {
		ostringstream msg;
		msg << "chol() - Error calling lapack::potrf(). Info = " << info;
		error(msg.str());
	}

	// Parse
	A_ = ublas::triangular_adaptor<CM, ublas::upper>(A_);
	return A_;
}

/// LU factorization.
/// Call PLU = lu( A ). Returns a permutation matrix P, a lower triangular
/// (trapezoidal) matrix L and an upper-triangular (trapezoidal) matrix X
/// such that A = P * L * U.
template<class E, class T = typename E::value_type>
tuple< imat, Mat<T>, Mat<T> > lu( const matrix_interface<E>& A )
{
	typedef ublas::matrix<T,ublas::column_major> CM;
	CM A_ = A;

	// Call parameters
	int M     = A_.size1();
	int N     = A_.size2();
	int lda   = max( 1, M );
	int mdim  = min(M,N);
	ivec ipiv(mdim);
	int info  = lapack::getrf( M, N, &A_(0,0), lda, &ipiv(0) );
	if( info < 0 ) {
		ostringstream msg;
		msg << "lu() - Error calling lapack::getrf(). Info = " << info;
		error(msg.str());
	}

	// Parse L
	Mat<T> L = ublas::triangular_adaptor<CM, ublas::lower>( A_ );
	L        = L( range(0,L.size1()), range(0,mdim) );

	// Parse U
	diag(L)  = ones( mdim );
	Mat<T> U = ublas::triangular_adaptor<CM, ublas::upper>( A_ );
	U        = U( range(0,mdim), range(0,U.size2()) );

	// Parse P
	imat P = eye( M );
	for( int i = 0; i < mdim; i++ ) {
		imat Pi = eye(M);			// i-th permutation matrix
		int j   = ipiv(i)-1;		// index with which the i-th row permutes
		Pi(i,i) = 0;
		Pi(i,j) = 1;
		Pi(j,j) = 0;
		Pi(j,i) = 1;
		P = P * Pi;
	}

	// Return a tuple
	return tuple<imat,Mat<T>,Mat<T>>{ P, L, U };
}
/// @}

//------------------------------------------------------------------------------
/// @name Inverses

/// Inverse of a matrix
template<class E, class T = typename E::value_type>
Mat<T> inv( const matrix_interface<E>& A )
{
	assert( A.size1() == A.size2() );
	int N = A.size1();
	Mat<T> I = eye(N);
	return linsolve( A, I );
}

/// Pseudo-inverse.
/// Call: Y = pinv( X, tol ).
template<class E, class T = typename E::value_type>
Mat<T> pinv( const matrix_interface<E>& A, double tol=-1 )
{
	if( tol == -1 )
		tol = eps * norm_op(A) * max(A.size1(),A.size2());
	auto USV  = svd( A );
	Mat<T> U  = get<0>(USV);
	rmat S    = get<1>(USV);
	Mat<T> V  = get<2>(USV);
	int N     = min( A.size1(), A.size2() );
	rmat Sp(A.size2(),A.size1());
	for( int n = 0; n < N; n++ ) {
		double Snn = S(n,n);
		if( Snn > tol )
			Sp(n,n) = 1 / S(n,n);
		else
			Sp(n,n) = 0;
	}
	return noproxy(V * Sp) * adjoint(U);
}
/// @}

//------------------------------------------------------------------------------
/// @name Functionals

/// Trace
template<class E, class T = typename E::value_type>
T trace( const matrix_interface<E>& A )
{
	return sum(diag(A));
}

/// Determinant
template<class E, class T = typename E::value_type>
T det( const matrix_interface<E>& A )
{
	// Assert dimensions
	assert( A.size1() == A.size2() );

	// LU factorization
	auto PLU = lu( A );
	imat P   = get<0>( PLU );
	Mat<T> U = get<2>( PLU );

	// Determinant modulus
	T mod    = prod( diag(U) );

	// Determinant sign
	int sign = 1;
	int I    = P.size1();
	for( int i = 0; i < I; i++ ) {
		int j = max_index( P.column(i) );		// Index with which the i-th row was exchanged
		if( j != i ) {
			ivec rtmp = P.row(j);
			P.row(j)  = P.row(i);
			P.row(i)  = rtmp;
			sign *= -1;
		}
	}

	// Return
	return sign * mod;
}

/// Operator norm
template<class E>
double norm_op( const matrix_interface<E> &A )
{
	typedef typename E::value_type T;
	auto USV = svd( A );
	Mat<T> U = get<0>(USV);
	rmat S   = get<1>(USV);
	Mat<T> V = get<2>(USV);
	return max( diag(S) );
}

/// Operator norm
template<class E>
double cond( const matrix_interface<E> &A )
{
	typedef typename E::value_type T;
	auto USV = svd( A );
	Mat<T> U = get<0>(USV);
	rmat S   = get<1>(USV);
	Mat<T> V = get<2>(USV);
	return max( diag(S) ) / min( diag(S) );
}
/// @}

//------------------------------------------------------------------------------
/// @name Functions on square matrices

/// Apply a scalar function to a matrix as a matrix function
template< class F, class E>
cmat matrix_function( const F& fun, const matrix_interface<E>& A )
{
	// Diagonal decomposition
	auto dv  = eig( A );
	cmat d   = get<0>( dv );
	cmat V   = get<1>( dv );

	// Apply the function
	int N = d.size1();
	for( int n = 0; n < N; n++ )
		d(n,n) = fun(d(n,n));
	cmat R = trans( linsolve( trans(V), trans(V*d) ) );

	return R;
}

/// Converts a scalar function 'fun' into a matrix one 'mfun'
#define _MATRIX_FUNCTION( fun )													\
template<class E>																\
cmat m##fun( const matrix_interface<E>& A )										\
{																				\
	auto functor = []( complex x ){ return fun(x); };							\
	return matrix_function( functor, A );										\
}

/// Matrix exponential
_MATRIX_FUNCTION( exp )

/// Matrix logarithm
_MATRIX_FUNCTION( log )

/// Matrix square root
_MATRIX_FUNCTION( sqrt )

/// Matrix power
template<class E, class T = typename E::value_type>
cmat mpow( const matrix_interface<E>& A, double p )
{
	// Init
	assert( A.size1() == A.size2() );
	int I = A.size1();
	cmat B;

	// Integer power
	if( mod(p,1) == 0 ) {
		B = eye(I);
		if( p >= 0 )
			for( int i = 0; i < p; i++ )
				B = B*A;
		else {
			Mat<T> Ai = inv( A );
			for( int i = 0; i < -p; i++ )
				B = B*Ai;
		}
	}
	// Real power
	else {
		auto fun = [p]( complex x ){ return pow(x,p); };
		B = matrix_function( fun, A );
	}
	return B;
}
/// @}

/// @}
}

#endif
