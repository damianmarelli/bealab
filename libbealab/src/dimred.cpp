#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOGCLASSES

#include <bealab/extensions/other/dimred.hpp>
#include <GClasses/GManifold.h>

using namespace GClasses;

namespace bealab
{
namespace dimred
{
//------------------------------------------------------------------------------
// Function multidimensional_scalling()
//------------------------------------------------------------------------------
Vec<rvec> multidimensional_scalling( const rmat& distances, int target_dimensions )
{
	assert( distances.size1() == distances.size2() );

	// Input matrix
	int I = distances.size1();
	GMatrix A(I,I);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < I; j++ )
			A[i][j] = distances(i,j);

	// Call the method
	GRand prand(0);
	GMatrix* pB = GManifold::multiDimensionalScaling( &A, target_dimensions,
			&prand, false );
	Holder<GMatrix> hB(pB);
	GMatrix& B = *pB;

	// Output data
	Vec<rvec> rv(I);
	for( int i = 0; i < I; i++ ) {
		rv(i).resize(target_dimensions);
		for( int j = 0; j < target_dimensions; j++ )
			rv(i)(j) = B[i][j];
	}
	return rv;
}

//------------------------------------------------------------------------------
// Class manifold_learning_b
//------------------------------------------------------------------------------
manifold_learning_b::manifold_learning_b( int nn, int td ) :
	n_neighbors(nn), target_dims(td)
{
	prng    = new GRand(0);
//		pmethod = new ml_method( nn, td, reinterpret_cast<GRand*>(prng) );
}

manifold_learning_b::~manifold_learning_b()
{
	delete reinterpret_cast<GRand*>(prng);
//		delete pmethod;
}

Vec<rvec> manifold_learning_b::run( const Vec<rvec>& points )
{
	// Input data
	int I = points.size();
	int J = points(0).size();
	GMatrix A(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			A[i][j] = points(i)(j);

	// Run manifold learning method
	GMatrix* pB = reinterpret_cast<GMatrix*>( run_mlearning( &A ) );
	Holder<GMatrix> hB(pB);
	GMatrix& B = *pB;

	// Output data
	Vec<rvec> rv(I);
	for( int i = 0; i < I; i++ ) {
		rv(i).resize(target_dims);
		for( int j = 0; j < target_dims; j++ )
			rv(i)(j) = B[i][j];
	}
	return rv;
}

//------------------------------------------------------------------------------
// Class manifold_sculpting
//------------------------------------------------------------------------------
manifold_sculpting::manifold_sculpting( int nn, int td ) :
	manifold_learning_b(nn,td)
{
	pmethod = new GManifoldSculpting( nn, td, reinterpret_cast<GRand*>(prng) );
}

manifold_sculpting::~manifold_sculpting()
{
	delete reinterpret_cast<GManifoldSculpting*>(pmethod);
}

void* manifold_sculpting::run_mlearning( const void* cvpA )
{
	void* vpA   = const_cast<void*>(cvpA);
	GMatrix* pA = reinterpret_cast<GMatrix*>(vpA);
	return reinterpret_cast<GManifoldSculpting*>(pmethod)->doit( *pA );
}

}
}
#endif
