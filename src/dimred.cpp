/*******************************************************************************
 * This software is licensed under the BSD 3-Clause License with the possibility to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
 * You may not use this work except in compliance with the License.
 * You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
 * Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the file License for the specific language governing permissions and limitations under the License. 
 * If you wish to obtain a commercial license, please contact the authors via e-mail.
 *
 * Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)
 *******************************************************************************/
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
// Function multidimensional_scaling() (GClasses version)
//------------------------------------------------------------------------------
vec<rvec> multidimensional_scaling_gc( const rmat& distances, int target_dimensions )
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
	vec<rvec> rv(I);
	for( int i = 0; i < I; i++ ) {
		rv(i).resize(target_dimensions);
		for( int j = 0; j < target_dimensions; j++ )
			rv(i)(j) = B[i][j];
	}
	return rv;
}

//------------------------------------------------------------------------------
// Class multidimensional_scaling
//------------------------------------------------------------------------------
multidimensional_scaling::multidimensional_scaling( const rmat& distances, int target_dimensions ) :
	D(distances), dim(target_dimensions) {}

vec<rvec> multidimensional_scaling::run()
{
	// Eigenvalue decomposition
	int N    = D.size1();
	rmat D2  = element_prod(D,D);
	rmat H   = eye(N) - ones(N,N)/(double)N;
	rmat B   = -1./2 * H * rmat(D2 * H);
	cmat M, E;
	tie(E,M) = eig( B );

	// Choose indexes
	rvec e     = real(diag(E));
	rvec e0    = element_prod( e, ( e > 0 ) );
	ivec sidxs = flip( sort_indices( e0 ) );
	ivec idxs  = sidxs( range(0,dim) );

	// Prepare output
	rmat X = real( M * diag(sqrt(e0)) );
	vec<rvec> rv(N);
	for( int n = 0; n < N; n++ )
		rv(n) = X.row(n)( indirect(idxs) );

	// Store eigenvalues
	eigenvalues = e;
	ordered_eigenvalues = e0(indirect(sidxs));

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

vec<rvec> manifold_learning_b::run( const vec<rvec>& points )
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
	vec<rvec> rv(I);
	for( int i = 0; i < I; i++ ) {
		rv(i).resize(target_dims);
		for( int j = 0; j < target_dims; j++ )
			rv(i)(j) = B[i][j];
	}
	return rv;
}

//------------------------------------------------------------------------------
// Class isomap
//------------------------------------------------------------------------------
isomap::isomap( int nn, int td ) :
	manifold_learning_b(nn,td)
{
	pmethod = new GIsomap( nn, td, reinterpret_cast<GRand*>(prng) );
}

isomap::~isomap()
{
	delete reinterpret_cast<GIsomap*>(pmethod);
}

void* isomap::run_mlearning( const void* cvpA )
{
	void* vpA   = const_cast<void*>(cvpA);
	GMatrix* pA = reinterpret_cast<GMatrix*>(vpA);
	return reinterpret_cast<GIsomap*>(pmethod)->doit( *pA );
}

//------------------------------------------------------------------------------
// Class lle
//------------------------------------------------------------------------------
lle::lle( int nn, int td ) :
	manifold_learning_b(nn,td)
{
	pmethod = new GLLE( nn, td, reinterpret_cast<GRand*>(prng) );
}

lle::~lle()
{
	delete reinterpret_cast<GLLE*>(pmethod);
}

void* lle::run_mlearning( const void* cvpA )
{
	void* vpA   = const_cast<void*>(cvpA);
	GMatrix* pA = reinterpret_cast<GMatrix*>(vpA);
	return reinterpret_cast<GLLE*>(pmethod)->doit( *pA );
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

//------------------------------------------------------------------------------
// Class breadth_first_unfolding
//------------------------------------------------------------------------------
breadth_first_unfolding::breadth_first_unfolding( int a, int nn, int td ) :
	manifold_learning_b(nn,td)
{
	pmethod = new GBreadthFirstUnfolding( a, nn, td, reinterpret_cast<GRand*>(prng) );
}

breadth_first_unfolding::~breadth_first_unfolding()
{
	delete reinterpret_cast<GBreadthFirstUnfolding*>(pmethod);
}

void* breadth_first_unfolding::run_mlearning( const void* cvpA )
{
	void* vpA   = const_cast<void*>(cvpA);
	GMatrix* pA = reinterpret_cast<GMatrix*>(vpA);
	return reinterpret_cast<GBreadthFirstUnfolding*>(pmethod)->doit( *pA );
}

}
}
#endif
