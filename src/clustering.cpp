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

namespace bealab
{
//---------------------------------------------------------------------------------------
// K-means clustering algorithm (alternative implementation):
//
// 	db      = database of training vectors
//  K       = number of centroids
//  passign = pointer to a vector indicating the centroid assigned to each training vector
//  palloc  = pointer to a vector indicating the number of training vector assigned to each centroid
//  pcount  = pointer to an integer counting the number of iterations
//---------------------------------------------------------------------------------------
vec<rvec> kmeans( const vec<rvec> &db, int K, ivec *passign, ivec *palloc, int* pcount )
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
