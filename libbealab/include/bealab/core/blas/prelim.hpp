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
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
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
