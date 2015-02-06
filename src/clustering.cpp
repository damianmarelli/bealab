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
#include <bealab/scilib/stats.hpp>

void clustr ( double x[], double d[], double dev[], int b[], double f[],
  int e[], int i, int j, int n, int nz, int k, int *pcount );

namespace bealab
{
//---------------------------------------------------------------------------------------
// K-means clustering algorithm
//
// 	db      = database of training vectors
//  K       = number of centroids
//  passign = pointer to a vector indicating the centroid assigned to each training vector
//  palloc  = pointer to a vector indicating the number of training vector assigned to each centroid
//  pcount  = pointer to an integer counting the number of iterations
//---------------------------------------------------------------------------------------
vec<rvec> kmeans( const vec<rvec> &pts, int K, ivec *passig, ivec *palloc, int *pcount )
{
	int N = pts.size();			// Points
	int M = pts(0).size();		// Dimensions

	//  Generate the data.
	double *points = new double[N*M];
	for( int n = 0; n < N; n++ )
		for( int m = 0; m < M; m++ )
			points[n+m*N] = pts(n)(m);

	//  Initialize the cluster centers randomly
	double *clusters = new double[K*M];
	for( int k = 0; k < K; k++ ) {
		int n = floor( N * randu() );
		for( int m = 0; m < M; m++ )
			clusters[k+m*K] = pts(n)(m);
	}

	//  Compute the clusters.
	double *deviation = new double[K];	// Return the squared deviation whithin each cluster
	int *assignments  = new int[N];		// Cluster of each point
	double *ws        = new double[N];	// Workspace ???
	int *alloc        = new int[K];		// Number of points on each cluster
	int Pmin          = 1;				// Minimum number of points per cluster
	int Kmax          = K;				// Maximum number of clusters
	clustr( points, clusters, deviation, assignments, ws, alloc, N, M, K, Pmin, Kmax, pcount );

	// Form the assignment vector
	if( passig != NULL ) {
		passig->resize( N );
		for( int n = 0; n < N; n++ )
			(*passig)(n) = assignments[n] - 1;
	}

	// Form the point allocation vector
	if( palloc != NULL ) {
		palloc->resize( K );
		for( int k = 0; k < K; k++ )
			(*palloc)(k) = alloc[k];
	}

	// Form the return vector
	vec<rvec> rv(K);
	for( int k = 0; k < K; k++ ) {
		rvec v(M);
		for( int m = 0; m < M; m++ )
			v(m) = clusters[k+m*K];
		rv(k) = v;
	}

	// Free
	delete[] points;
	delete[] clusters;
	delete[] deviation;
	delete[] assignments;
	delete[] ws;
	delete[] alloc;

	return rv;
}

//---------------------------------------------------------------------------------------
// K-means clustering algorithm (alternative implementation):
//
// 	db      = database of training vectors
//  K       = number of centroids
//  passign = pointer to a vector indicating the centroid assigned to each training vector
//  palloc  = pointer to a vector indicating the number of training vector assigned to each centroid
//  pcount  = pointer to an integer counting the number of iterations
//---------------------------------------------------------------------------------------
vec<rvec> kmeans1( const vec<rvec> &db, int K, ivec *passign, ivec *palloc, int* pcount )
{
	// Initial centroids
	int N = db.size();
	vec<rvec> centroids(K);
	for( int k = 0; k < K; k++ ) {
		centroids(k) = db( rand(uniform_discrete(0,N-1)) );
		for( int l = 0; l < k; l++ )
			if( sum( centroids(l) == centroids(k) ) ) {
				k--;
				continue;
			}
	}

	// Main loop
	ivec assign = -1*ones(N);
	ivec assign0;
	ivec alloc(K);
	int count = 0;
	do {

		// Assign each training vector to the nearest centroid
		assign0 = assign;
		for( int n = 0; n < N; n++ ) {
			rvec tvec      = db(n);
			double mindist = inf;
			int minidx     = -1;
			for( int k = 0; k < K; k++ ) {
				rvec diff = tvec-centroids(k);
				double dist = inner_prod( diff, diff );
				if( dist < mindist ) {
					mindist = dist;
					minidx  = k;
				}
			}
			assign(n) = minidx;
		}

		// Recompute the centroids
		alloc.clear();
		for( int k = 0; k < K; k++ )
			centroids(k).clear();
		for( int n = 0; n < N; n++ ) {
			int cidx = assign(n);
			alloc(cidx)++;
			centroids(cidx) += db(n);
		}
		for( int k = 0; k < K; k++ ) {
			if( alloc(k) == 0 )
				centroids(k) = db( rand(uniform_discrete(0,N-1)) );
			else
				centroids(k) /= alloc(k);
		}

		// Update iteration counter
		count++;

	} while( sum(abs(assign-assign0)) != 0 );

	// return values
	if( passign != NULL )
		*passign = assign;
	if( palloc != NULL )
		*palloc = alloc;
	if( pcount != NULL )
		*pcount = count;
	return centroids;
}

}
