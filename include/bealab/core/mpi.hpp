// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/core/mpi.hpp
/// Support for cluster parallelization.

#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOMPI

#ifndef _BEALAB_MPI_
#define	_BEALAB_MPI_

#include <bealab/core/blas/vector.hpp>
#include <boost/mpi.hpp>

namespace bealab
{
/// @defgroup mpi Message Passing Interface (MPI)
/// Support for cluster parallelization.
/// @{

// Use namespace
namespace mpi = boost::mpi;

/// Implements a for loop in parallel over a cluster.
/// It evaluates fun(i), for i=0,...,I-1, and return a vector with the results.
template<class T>
Vec<T> parallel_for( int I, const function<T(int)>& fun )
{
	// Parallel processing
	Vec<T> X(I);
	mpi::communicator com;
	#pragma omp parallel for
	for( int i = com.rank(); i < I; i += com.size() )
		X[i] = fun(i);

	// Collect the result
	vector<Vec<T>> XX;
	mpi::gather( com, X, XX, 0 );
	Vec<T> Y(I);
	if( com.rank() == 0 )
		for( int i = 0; i < I; i++ )
			Y[i] = XX[ i % com.size() ](i);

	// Share the result with the other nodes
	mpi::broadcast( com, Y, 0 );

	return Y;
}

/// @}
}
#endif
#endif
