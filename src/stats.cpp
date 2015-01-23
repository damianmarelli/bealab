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

// Random number generation ----------------------------------------------------
namespace _stats
{
// Random engine
std::mt19937 engine;
//std::random_device rd;
//std::default_random_engine engine;

// Normal distribution mapping for random number generation
std::normal_distribution<double> stats2random(
		const boost::math::normal_distribution<double>& dist )
{
	return std::normal_distribution<double>( dist.mean(), dist.standard_deviation() );
}

}

// Initializes the engine with a random seed (obtained from std::time()).
void randomize( unsigned int seed )
{
	_stats::engine.seed( seed );
}

// Multivariate normal distribution --------------------------------------------
// Probability density function
double pdf( const multivariate_normal& dist, const rvec& x )
{
	int N  = x.size();
	rvec m = dist.mean();
	rmat R = dist.covariance();
	assert( N == (int)m.size() );
	rvec e   = x-m;
	double y = inner_prod( e, linsolve( R, e ) );
	return 1/sqrt( pow(2*pi,N) * det(R) ) * exp( -1./2 * y );
}

#ifndef BEALAB_NOFORTRAN

// Cumulative distribution function over a square region
//double cdf( const multivariate_normal& dist, const rvec& a, const rvec& b )
//{
//	assert( a.size() == b.size() );
//	assert( a.size() == dist.size() );
//
//	if( a.size() == 1 ) {
//		normal sn( dist.mean()(0), dist.covariance()(0,0) );
//		return cdf( sn, b(0) ) - cdf( sn, a(0) );
//	}
//	else {
//		function<double(const rvec&)> p = [&dist]( const rvec& y ) { return pdf( dist, y ); };
//		return integral( p, a, b );
//	}
//}
extern "C" {
	void mvndst_( int* N, double* lower, double* upper, int* infin,
			double* correl, int* maxpts, double* abseps, double* releps,
			double* error, double* value, int* inform );
}
double cdf( const multivariate_normal& dist, rvec a, rvec b )
{
	// Assertions
	assert( a.size() == dist.size() );
	assert( a.size() == b.size() );

	// Integration limits flags
	a -= dist.mean();
	b -= dist.mean();
	int N = dist.size();
	ivec infin = zeros(N);
	for( int n = 0; n < N; n++ )
		if( isinf(a(n)) )
			if( isinf(b(n)) )
				infin(n) = -1;
			else
				infin(n) = 0;
		else
			if( isinf(b(n)) )
				infin(n) = 1;
			else
				infin(n) = 2;

	// Correlation matrix
	const rmat& C = dist.covariance();
	double correl[N*(N-1)/2];
	for( int i = 0; i < N; i++ )
		for( int j = 0; j < i; j++ )
			correl[ j + ((i-1)*i)/2 ] = C(i,j);

	// Other call parameters
	int maxpts    = 1000*N;
	double abseps = N < 4 ? 1e-8 : 1e-4;
	double releps = 1;
	double err;
	double value;
	int inform;

	// Call Fortran routine
	mvndst_( &N, &a(0), &b(0), &infin(0), correl,
			&maxpts, &abseps, &releps, &err, &value, &inform );
	if( inform )
		error("Multivariate normal CDF - Error calling mvndst()");

	return value;
}

// Cumulative distribution function over a square region
double cdf( const multivariate_normal& dist, const rvec& b )
{
	assert( b.size() == dist.size() );
	rvec a = -inf * ones( b.size() );
	return cdf( dist, a, b );
}

#endif

}
