/// @file bealab/scilib/clustering.hpp
/// Clustering routines.

#ifndef _BEALAB_CLUSTERING_
#define	_BEALAB_CLUSTERING_

#include <bealab/core/blas.hpp>

namespace bealab
{
/// @defgroup clustering Clustering
/// Clustering routines.
/// @{

/// K-means algorithm over a vector of vector samples.
vec<rvec> kmeans( const vec<rvec>&, int, ivec* =NULL, ivec* =NULL, int* =NULL );
vec<rvec> kmeans1( const vec<rvec>&, int, ivec* =NULL, ivec* =NULL, int* =NULL );

/// @}
}

#endif
