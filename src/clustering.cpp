// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

#include <bealab/scilib/stats.hpp>

namespace bealab
{

/*
  Clustering function used by kmeans

    Given a matrix of I observations on J variables, the
    observations are allocated to N clusters in such a way that the
    within-cluster sum of squares is minimised.

  Parameters:

    Input, double X[I*J], the observed data.

    Input/output, double D[K*J], the cluster centers.
    On input, the user has chosen these.  On output, they have been
    updated.

    Output, double DEV[K], the sums of squared deviations
    of observations from their cluster centers.

    Output, int B[I], indicates the cluster to which
    each observation has been assigned.

    Workspace, double F[I].

    Output, int E[K], the number of observations assigned
    to each cluster.

    Input, int I, the number of observations.

    Input, int J, the number of variables.

    Input, int N, the number of clusters.

    Input, int NZ, the minimum number of observations
    which any cluster is allowed to have.

    Input, int K, the maximum number of clusters.
*/
static
void clustr ( double x[], double d[], double dev[], int b[], double f[],
  int e[], int i, int j, int n, int nz, int k, int *pcount )
{
  double big = 1.0E+10;
  double da;
  double db;
  double dc;
  double de;
  double fl;
  double fm;
  double fq;
  int ia;
  int ic;
  int id;
  int ie;
  int ig;
  int ih;
  int ii;
  int ij;
  int ik;
  int il;
  int in;
  int ip;
  int ir;
  int is;
  int it;
  int iu;
  int iw;
  int ix;
  int iy;

  for ( ia = 1; ia <= n; ia++ )
  {
    e[ia-1] = 0;
  }
/*
  For each observation, calculate the distance from each cluster
  center, and assign to the nearest.
*/
  for ( ic = 1; ic <= i; ic++ )
  {
    f[ic-1] = 0.0;
    da = big;

    for ( id = 1; id <= n; id++ )
    {
      db = 0.0;
      for ( ie = 1; ie <= j; ie++ )
      {
        dc = x[ic-1+(ie-1)*i] - d[id-1+(ie-1)*k];
        db = db + dc * dc;
      }

      if ( db < da )
      {
        da = db;
        b[ic-1] = id;
      }
    }
    ig = b[ic-1];
    e[ig-1] = e[ig-1] + 1;
  }
/*
  Calculate the mean and sum of squares for each cluster.
*/
  for ( ix = 1; ix <= n; ix++ )
  {
    dev[ix-1] = 0.0;
    for ( iy = 1; iy <= j; iy++ )
    {
      d[ix-1+(iy-1)*k] = 0.0;
    }
  }

  for ( ic = 1; ic <= i; ic++ )
  {
    ig = b[ic-1];
    for ( ih = 1; ih <= j; ih++ )
    {
      d[ig-1+(ih-1)*k] = d[ig-1+(ih-1)*k] + x[ic-1+(ih-1)*i];
    }
  }

  for ( ij = 1; ij <= j; ij++ )
  {
    for ( ii = 1; ii <= n; ii++ )
    {
      d[ii-1+(ij-1)*k] = d[ii-1+(ij-1)*k] / ( double ) e[ii-1];
    }
  }

  for ( ij = 1; ij <= j; ij++ )
  {
    for ( ik = 1; ik <= i; ik++ )
    {
      il = b[ik-1];
      da = x[ik-1+(ij-1)*i] - d[il-1+(ij-1)*k];
      db = da * da;
      f[ik-1] = f[ik-1] + db;
      dev[il-1] = dev[il-1] + db;
    }
  }

  for ( ik = 1; ik <= i; ik++ )
  {
    il = b[ik-1];
    fl = e[il-1];
    if ( 2.0 <= fl )
    {
      f[ik-1] = f[ik-1] * fl / ( fl - 1.0 );
    }
  }
/*
  Examine each observation in turn to see if it should be
  reassigned to a different cluster.
*/
  if(pcount != NULL) *pcount = 1;

  for ( ; ; )
  {
    iw = 0;

    for ( ik = 1; ik <= i; ik++ )
    {
      il = b[ik-1];
      ir = il;
/*
  If the number of cluster points is less than or equal to the
  specified minimum, NZ, then bypass this iteration.
*/
      if ( nz < e[il-1] )
      {
        fl = e[il-1];
        dc = f[ik-1];

        for ( in = 1; in <= n; in++ )
        {
          if ( in != il )
          {
            fm = e[in-1];
            fm = fm / ( fm + 1.0 );

            de = 0.0;
            for ( ip = 1; ip <= j; ip++ )
            {
              da = x[ik-1+(ip-1)*i] - d[in-1+(ip-1)*k];
              de = de + da * da * fm;
            }

            if ( de < dc )
            {
              dc = de;
              ir = in;
            }
          }
        }
/*
  Reassignment is made here if necessary.
*/
        if ( ir != il )
        {
          fq = e[ir-1];
          dev[il-1] = dev[il-1] - f[ik-1];
          dev[ir-1] = dev[ir-1] + dc;
          e[ir-1] = e[ir-1] + 1;
          e[il-1] = e[il-1] - 1;

          for ( is = 1; is <= j; is++ )
          {
            d[il-1+(is-1)*k] = ( d[il-1+(is-1)*k]
                             * fl - x[ik-1+(is-1)*i] ) / ( fl - 1.0 );
            d[ir-1+(is-1)*k] = ( d[ir-1+(is-1)*k]
                             * fq + x[ik-1+(is-1)*i] ) / ( fq + 1.0 );
          }

          b[ik-1] = ir;

          for ( it = 1; it <= i; it++ )
          {
            ij = b[it-1];

            if ( ij == il || ij == ir )
            {
              f[it-1] = 0.0;
              for ( iu = 1; iu <= j; iu++ )
              {
                da = x[it-1+(iu-1)*i] - d[ij-1+(iu-1)*k];
                f[it-1] = f[it-1] + da * da;
              }
              fl = e[ij-1];
              f[it-1] = f[it-1] * fl / ( fl - 1.0 );
            }
          }
          iw = iw + 1;
        }
      }
    }

    if(pcount != NULL) (*pcount)++;

/*
  If any reassignments were made on this pass, then do another pass.
*/
    if ( iw == 0 )
    {
      break;
    }
  }
  return;
}

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
	for( int k = 0; k < K; k++ )
		centroids(k) = db( rand(uniform_discrete(0,N-1)) );

	// Main loop
	ivec assign(N);
	ivec assign0;
	ivec alloc(K);
	int count = 0;
	do {

		// Assign each training vector to the nearest centroid
		assign0 = assign;
		for( int n = 0; n < N; n++ ) {
			rvec tvec   = db(n);
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

	} while( sum(assign-assign0) != 0 );

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
