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

/// Enable a template definition if _expression::value_type == _type
#define ENABLE_IF_VALUE_TYPE( _expression_, _type_ ) 							\
	typename enable_if< 														\
		is_same< typename _expression_::value_type, _type_ >::value				\
	>::type* = 0

/// @name Convert expression templates into temporaries
template<class> class vector_interface;
template<class> class matrix_interface;
template<class value_type> class vectorx;
template<class value_type> class matrixx;

template<class T>
T temporary( const T& x )
{
	return x;
};

template<class E>
vector_interface<vectorx<typename E::value_type>>
temporary( const vector_interface<E>& x )
{
	return x;
}

template<class E>
matrix_interface<matrixx<typename E::value_type>>
temporary( const matrix_interface<E>& x )
{
	return x;
}
/// @}

/// @}
}
#endif
