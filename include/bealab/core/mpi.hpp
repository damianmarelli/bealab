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
