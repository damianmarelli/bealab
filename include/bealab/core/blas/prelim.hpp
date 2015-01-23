// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/core/blas/prelim.hpp
/// BLAS preliminaries

#ifndef _BEALAB_BLAS_PRELIM_
#define	_BEALAB_BLAS_PRELIM_

// uBlas config
#ifdef NDEBUG
#define BOOST_UBLAS_NDEBUG
#endif

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <bealab/core/prelim.hpp>

namespace bealab
{
/// @defgroup blas_prelim Preliminaries
/// BLAS preliminaries
/// @{

// Import uBLAS members
namespace ublas = boost::numeric::ublas;
using ublas::slice;
using ublas::range;
using ublas::lower;
using ublas::upper;

/// Indirect indexing
class indirect : public ublas::indirect_array<> {
public:
	indirect() {}
	indirect( const ublas::vector<long unsigned int> &v ) :
		indirect_array<>( v.size(), v.data() ) {}
};

/// Mask indexing
class mask : public indirect {
public:
	mask( const ublas::vector<bool> &msk )
	{
		int N = msk.size();
		int M = ublas::sum( ublas::vector<int>(msk) );
		ublas::vector<long unsigned int> idxs(M);
		for( int n = 0, m = 0; n < N; n++ )
			if( msk(n) == true )
				idxs(m++) = n;
		*dynamic_cast<indirect*>(this) = indirect(idxs);
	}
};

/// Default conversion of an expression template into a temporary
template<class T>
T noproxy( const T& x )
{
	return x;
};

/// Check for scalar types (i.e., bool, int, double complex)
template< class T >
struct is_scalar : std::integral_constant<
	   bool,
	   std::is_same<bool, typename std::remove_cv<T>::type>::value  ||
	   std::is_same<int, typename std::remove_cv<T>::type>::value  ||
	   std::is_same<double, typename std::remove_cv<T>::type>::value  ||
	   std::is_same<complex, typename std::remove_cv<T>::type>::value
   > {};

/// @}
}
#endif
