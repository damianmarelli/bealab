#include <bealab/scilib/linalg/lapack.hpp>

//------------------------------------------------------------------------------
// Lapack Fortran declarations
//------------------------------------------------------------------------------
namespace bealab { namespace lapack {

extern "C" {

	/// Called to choose problem-dependent parameters for the local environment
	int ilaenv_( int* ispec, const char* name, const char* opts, int* n1, int* n2, int* n3, int* n4 );

	/// Singular value decomposition (real)
	void dgesdd_(char const* jobz, int* M, int* N,
			double* A, int* lda,
			double* s, double* U, int* ldu,
			double* VT, int* ldvt,
			double* work, int* lwork, int* iwork, int* info );

	/// Singular value decomposition (complex)
	void zgesdd_(char* jobz, int* M, int* N,
			complex* A, int* lda,
			double* s, complex* U, int* ldu,
			complex* VT, int* ldvt,
			complex* work, int* lwork, double* rwork, int* iwork, int* info );

	/// Least-squares solution (real)
	void dgelsd_( int* M, int* N, int* Nrhs, double* A, int* lda, double* B,
			int* ldb, double* s, double* rcond, int* rank, double* work,
			int* lwork, int* iwork, int* info );

	/// Least-squares solution (complex)
	void zgelsd_( int* M, int* N, int* Nrhs, complex* A, int* lda, complex* B,
			int* ldb, double* s, double* rcond, int* rank, complex* work,
			int* lwork, double* rwork, int* iwork, int* info );

	/// Solution to a linear system of equations (real)
	void dgesv_( int* N, int* Nrhs, double* A, int* lda, int* ipiv, double* B,
			int* ldb, int* info );

	/// Solution to a linear system of equations (real)
	void zgesv_( int* N, int* Nrhs, complex* A, int* lda, int* ipiv, complex* B,
			int* ldb, int* info );

	/// Eigenvalue decomposition (real)
	void dgeev_( char* jobl, char* jobr, int* N, double* A, int* lda, double* WR,
			double* WI, double* VL, int* ldvl, double* VR, int* ldvr, double* work,
			int* lwork, int* info );

	/// Eigenvalue decomposition (imag)
	void zgeev_( char* jobl, char* jobr, int* N, complex* A, int* lda, complex* W,
			complex* VL, int* ldvl, complex* VR, int* ldvr, complex* work,
			int* lwork, double* rwork, int* info );

	/// Cholesky factorization (real)
	void dpotrf_( char* uplo, int* N, double* A, int* lda, int* info );

	/// Cholesky factorization (imag)
	void zpotrf_( char* uplo, int* N, complex* A, int* lda, int* info );

	/// LU factorization (real)
	void dgetrf_( int* M, int* N, double* A, int* lda, int* ipiv, int* info );

	/// LU factorization (imag)
	void zgetrf_( int* M, int* N, complex* A, int* lda, int* ipiv, int* info );
}

//------------------------------------------------------------------------------
// Lapack C++ interface
//------------------------------------------------------------------------------

/// Called to choose problem-dependent parameters for the local environment
template<>
int ilaenv<double>( int ispec, string name, string opts, int n1, int n2, int n3, int n4 )
{
	name = "D" + name;
	return ilaenv_( &ispec, name.data(), opts.data(), &n1, &n2, &n3, &n4 );
}

template<>
int ilaenv<complex>( int ispec, string name, string opts, int n1, int n2, int n3, int n4 )
{
	name = "Z" + name;
	return ilaenv_( &ispec, name.data(), opts.data(), &n1, &n2, &n3, &n4 );
}

/// Singular value decomposition (real)
template<>
int gesdd( char jobz, int M, int N, double* A, int lda,
        double* s, double* U, int ldu, double* VT, int ldvt,
        double* work, int lwork, int* iwork )
{
	int info;
	dgesdd_( &jobz, &M, &N, A, &lda, s, U, &ldu, VT, &ldvt, work, &lwork, iwork, &info );
	return info;
}

/// Singular value decomposition (complex)
template<>
int gesdd( char jobz, int M, int N, complex* A, int lda,
        double* s, complex* U, int ldu, complex* VT, int ldvt,
        complex* work, int lwork, int* iwork )
{
	int info;
	int lrwork = 5*min(M,N)*min(M,N) + 7*min(M,N);
	double rwork[max(1,lrwork)];
	zgesdd_( &jobz, &M, &N, A, &lda, s, U, &ldu, VT, &ldvt, work, &lwork, rwork, iwork, &info );
	return info;
}

/// Least-squares solution (real)
template<>
int gelsd( int M, int N, int Nrhs, double* A, int lda,
		double* B, int ldb, double* s, double rcond, int& rank,
		double* work, int lwork, int* iwork )
{
	int info;
	dgelsd_( &M, &N, &Nrhs, A, &lda, B, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info );
	return info;

}

/// Least-squares solution (complex)
template<>
int gelsd( int M, int N, int Nrhs, complex* A, int lda,
		complex* B, int ldb, double* s, double rcond, int& rank,
		complex* work, int lwork, int* iwork )
{
	int info;
	int lrwork;
	int smlsiz = lapack::ilaenv<complex>( 9, "GELSD", "", 0, 0, 0, 0 );
	int nlvl   = max( 0, int( log2( min( M, N ) / (smlsiz+1) ) ) + 1 );
	if( M >= N )
		lrwork = 10*N + 2*N*smlsiz + 8*N*nlvl + 3*smlsiz*Nrhs + pow( smlsiz+1, 2 );
	else
		lrwork = 10*M + 2*M*smlsiz + 8*M*nlvl + 3*smlsiz*Nrhs + pow( smlsiz+1, 2 );
	double rwork[max(1,lrwork)];
	zgelsd_( &M, &N, &Nrhs, A, &lda, B, &ldb, s, &rcond, &rank, work, &lwork, rwork, iwork, &info );
	return info;

}

/// Solution to a linear system of equations (real)
template<>
int gesv( int N, int Nrhs, double* A, int lda, int* ipiv, double* B,
		int ldb )
{
	int info;
	dgesv_( &N, &Nrhs, A, &lda, ipiv, B, &ldb, &info );
	return info;
}

/// Solution to a linear system of equations (real)
template<>
int gesv( int N, int Nrhs, complex* A, int lda, int* ipiv, complex* B,
		int ldb )
{
	int info;
	zgesv_( &N, &Nrhs, A, &lda, ipiv, B, &ldb, &info );
	return info;
}

/// Eigenvalue decomposition (real) XXX Doesn't work
template<>
int geev( char jobl, char jobr, int N, double* A, int lda, complex* W, complex* VL,
		int ldvl, complex* VR, int ldvr, double* work, int lwork )
{
	// Lapack call
	int info;
	double WR[N], WI[N];
	double VL_[ldvl], VR_[ldvr];
	dgeev_( &jobl, &jobr, &N, A, &lda, WR, WI, VL_, &ldvl, VR_, &ldvr, work, &lwork, &info );

	if( lwork == -1 )
		return info;

	for( int n = 0; n < N; n++ )
		cout << "n: " << n << ", Real: " << WR[n] << ", Imag: " << WI[n] << endl;

	// Parse the results
	for( int n = 0; n < N; n++ ) {

		// Eigenvalues
		if( WI[n] == 0 ) {
			W[n]   = WR[n];
		}
		else {
			W[n]   = WR[n] + i * WI[n];
			W[n+1] = WR[n] - i * WI[n];
		}

		// Left eigenvectors
		if( jobl == 'V' ) {

			// If real eigenvalue: copy n-th column of VL_ to the n-th column of VL
			if( WI[n] == 0 ) {
				for( int m = 0; m < N; m++ )
					VL[m+n*N] = VL_[m+n*N];
			}

			// If complex eigenvalue: VL(:,n) = VL(:,n) + i*VL(:,n+1) and
			//                        VL(:,n) = VL(:,n) - i*VL(:,n+1)
			else {
				for( int m = 0; m < N; m++ ) {
					VL[m+n*N]     = VL_[m+n*N] + i * VL_[m+(n+1)*N];
					VL[m+(n+1)*N] = VL_[m+n*N] - i * VL_[m+(n+1)*N];
				}

			}
		}

		// Right eigenvectors
		if( jobr == 'V' ) {

			// If real eigenvalue: copy n-th column of VR_ to the n-th column of VR
			if( WI[n] == 0 ) {
				for( int m = 0; m < N; m++ )
					VR[m+n*N] = VR_[m+n*N];
			}

			// If complex eigenvalue: VR(:,n) = VR(:,n) + i*VR(:,n+1) and
			//                        VR(:,n) = VR(:,n) - i*VR(:,n+1)
			else {
				for( int m = 0; m < N; m++ ) {
					VR[m+n*N]     = VR_[m+n*N] + i * VR_[m+(n+1)*N];
					VR[m+(n+1)*N] = VR_[m+n*N] - i * VR_[m+(n+1)*N];
				}

			}
		}

		// If this eigenvalue was complex, skip the next
		if( WI[n] != 0 )
			n++;
	}

	return info;
}

/// Eigenvalue decomposition (complex)
template<>
int geev( char jobl, char jobr, int N, complex* A, int lda, complex* W, complex* VL,
		int ldvl, complex* VR, int ldvr, complex* work, int lwork )
{
	int info;
	double rwork[2*N];
	zgeev_( &jobl, &jobr, &N, A, &lda, W, VL, &ldvl, VR, &ldvr, work, &lwork,
			rwork, &info );
	return info;
}

/// Cholesky factorization (real)
template<>
int potrf( char uplo, int N, double* A, int lda )
{
	int info;
	dpotrf_( &uplo, &N, A, &lda, &info );
	return info;
}

/// Cholesky factorization (complex)
template<>
int potrf( char uplo, int N, complex* A, int lda )
{
	int info;
	zpotrf_( &uplo, &N, A, &lda, &info );
	return info;
}

template<>
int getrf( int M, int N, double* A, int lda, int* ipiv )
{
	int info;
	dgetrf_( &M, &N, A, &lda, ipiv, &info );
	return info;
}

template<>
int getrf( int M, int N, complex* A, int lda, int* ipiv )
{
	int info;
	zgetrf_( &M, &N, A, &lda, ipiv, &info );
	return info;
}

}}
