// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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

/// @}
}

#endif
