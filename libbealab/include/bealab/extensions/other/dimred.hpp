/// @file bealab/extensions/other/dimred.hpp
/// Dimensionality reduction methods.

#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOGCLASSES

#ifndef _BEALAB_DIMRED_
#define _BEALAB_DIMRED_

#include <bealab/scilib.hpp>

namespace bealab
{
namespace dimred
{
/// @defgroup dimred Dimensionality reduction
/// Dimensionality reduction methods.
/// @{

/// Compute the N principal components of the vectors in V
template<class T=double>
class principal_component_analysis {

	int N;

public:

	Vec<Vec<T>> principal_components;
	Vec<T> approximation;
	Vec<T> coefficients;

	principal_component_analysis( const Vec<Vec<T>>& x, int N_ ) :
		N(N_), principal_components(N_)
	{
		// Convert to a matrix (X)
		int I = x(0).size();
		int J = x.size();
		Mat<T> X(I,J);
		for( int j = 0; j < J; j++ )
			X.column(j) = x(j);

		// Compute principal components of X
		Mat<T> XX = X * trans(X);
		auto DU   = eig( XX );
		rvec d    = real(diag( get<0>(DU) ));
		rmat U    = real(get<1>(DU));
		Mat<T> PCS(I,N);
		for( int n = 0; n < N; n++ )
			PCS.column(n) = U.column(n);

		// Convert to an array of vectors (pcs)
		for( int n = 0; n < N; n++ )
			principal_components(n) = PCS.column(n);
	}

	rvec approximate( const Vec<T>& x )
	{
		// Compute the coefficients
		coefficients.resize(N);
		for( int n = 0; n < N; n++ )
			coefficients(n) = inner_prod( x, principal_components(n) );

		// Compute the best approximation
		approximation = zeros(x.size());
		for( int n = 0; n < N; n++ )
			approximation += coefficients(n) * principal_components(n);

		return approximation;
	}
};

/// Base for all manifold learning methods
class manifold_learning_b {
protected:

	void* pmethod = 0;
	void* prng;

	virtual
	void* run_mlearning( const void* pA ) = 0;

public:

	int n_neighbors;
	int target_dims;

	/// Constructor
	manifold_learning_b( int nn, int td );

	/// Destructor
	virtual
	~manifold_learning_b();

	/// Run the manifold pearning algorithm
	Vec<rvec> run( const Vec<rvec>& points );
};

///// Isomap
//struct isomap : public manifold_learning_b<GIsomap>
//{
//	using manifold_learning_b::manifold_learning_b;
//};
//
///// Locally linear embedding
//struct lle : public manifold_learning_b<GLLE>
//{
//	using manifold_learning_b::manifold_learning_b;
//};

/// Manifold sculpting
//struct manifold_sculpting : public manifold_learning_b<GManifoldSculpting>
struct manifold_sculpting : public manifold_learning_b
{
	manifold_sculpting( int nn, int td );
	~manifold_sculpting();
	void* run_mlearning( const void* cvpA );
};

///// Breadth first unfolding
//struct breadth_first_unfolding : public manifold_learning_b<GBreadthFirstUnfolding>
//{
//	breadth_first_unfolding( int a, int nn, int td )
//	{
//		n_neighbors = nn;
//		target_dims = td;
//		prng        = new GRand(0);
//		pmethod     = new GBreadthFirstUnfolding( a, nn, td, reinterpret_cast<GRand*>(prng) );
//	}
//};

/// Multi-dimensional scalling method
Vec<rvec> multidimensional_scalling( const rmat& distances, int target_dimensions );

/// @}
}
}
#endif
#endif
