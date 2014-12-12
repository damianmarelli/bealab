/// @file bealab/scilib/linalg/lapack.hpp
/// Interfaces

#ifndef _BEALAB_LAPACK_
#define	_BEALAB_LAPACK_

#include <bealab/core/prelim.hpp>

namespace bealab
{
namespace lapack
{

/// @defgroup linalg_lapack Lapack interface
/// Interfaces to the LAPACK libraries
/// @{

/// Called to choose problem-dependent parameters for the local environment
template<class T>
int ilaenv( int ispec, string name, string opts, int n1, int n2, int n3, int n4 );

/// Singular value decomposition
template<class T>
int gesdd( char jobz, int M, int N, T* A, int lda, double* s, T* U, int ldu,
		T* VT, int ldvt, T* work, int lwork, int* iwork );

/// Least-squares solution
template<class T>
int gelsd( int M, int N, int Nrhs, T* A, int lda, T* B, int ldb, double* s,
		double rcond, int& rank, T* work, int lwork, int* iwork );

/// Solution to a linear system of equations
template<class T>
int gesv( int N, int Nrhs, T* A, int lda, int* ipiv, T* B, int ldb );

/// Eigenvalue decomposition
template<class T>
int geev( char jobl, char jobr, int N, T* A, int lda, complex* W, complex* VL,
		int ldvl, complex* VR, int ldvr, T* work, int lwork );

/// Cholesky factorization
template<class T>
int potrf( char uplo, int N, T* A, int lda );

/// LU factorization
template<class T>
int getrf( int M, int N, T* A, int lda, int* ipiv );

/// @}
}
}
#endif

