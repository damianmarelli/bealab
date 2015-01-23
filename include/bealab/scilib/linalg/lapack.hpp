// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/scilib/linalg/lapack.hpp
/// Interfaces

#ifndef _BEALAB_LAPACK_
#define	_BEALAB_LAPACK_

#include <bealab/core/prelim.hpp>

namespace bealab
{
/// Lapack interface
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

