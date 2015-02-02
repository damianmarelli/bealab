/// @file bealab/scilib/stats.hpp
/// Analysis of statistical distributions and random number generation.

#ifndef _BEALAB_STATS_
#define	_BEALAB_STATS_

#include <random>
#include <ctime>
#include <boost/math/distributions/normal.hpp>
#include <bealab/core/blas.hpp>
#include <bealab/scilib/calculus.hpp>
#include <bealab/scilib/linalg.hpp>

namespace bealab
{
/// @defgroup stats Statistics
/// Analysis of statistical distributions and random number generation.
/// @{

//------------------------------------------------------------------------------
/// @name Univariate distributions
typedef boost::math::normal_distribution<double> normal;						///< Normal distribution
typedef std::uniform_real_distribution<double> uniform;							///< Uniform continuous distribution
typedef std::uniform_int_distribution<int> uniform_discrete;					///< Uniform discrete distribution
typedef std::bernoulli_distribution bernoulli;									///< Bernoulli distribution
/// @}

//------------------------------------------------------------------------------
/// @name Multivariate normal distribution

/// Distribution class
class multivariate_normal {

	rvec _mean;																	///< Mean vector
	rmat _covariance;															///< Covariance matrix
	rmat _covariance_factor;													///< Covariance factor (i.e., X such that C = X * trans(X))

public:

	/// Constructor
	multivariate_normal( const rvec& m, const rmat& Cf ) : _mean(m), _covariance_factor(Cf)
	{
		assert( m.size()  == Cf.size1() );
	}

	/// Get the mean vector
	const rvec& mean() const { return _mean; }

	/// Get the covariance matrix
	const rmat covariance() const
	{
		if( _covariance.size1() > 0 )
			return _covariance;
		else
			return _covariance_factor * trans(_covariance_factor);
	}

	/// Set the covariance matrix
	void covariance( const rmat& cov )
	{
		_covariance = cov;
	}

	/// Dimension of the random vector
	uint size() const { return _mean.size(); }

	/// Generate a random vector
	template<class generator>
	rvec operator()( generator& gen )
	{
		int M = size();
		rvec rv(M);
		std::normal_distribution<double> dist( 0, 1 );
		for( int m = 0; m < M; m++ )
			rv(m) = dist(gen);
		return _covariance_factor * rv + _mean;
	}
};

/// Probability density function
double pdf( const multivariate_normal& dist, const rvec& x );

/// Cumulative distribution function over a square region
double cdf( const multivariate_normal& dist, rvec a, rvec b );

/// Cumulative distribution function over a square region
double cdf( const multivariate_normal& dist, const rvec& b );

/// @}

//------------------------------------------------------------------------------
/// @name Random number generation

/// Private namespace for internal use.
namespace _stats
{

/// Random engine
extern std::mt19937 engine;
//extern std::default_random_engine engine;

/// Normal distribution mapping for random number generation
std::normal_distribution<double> stats2random(
		const boost::math::normal_distribution<double>& dist );

/// multivariate normal distribution mapping for random number generation
inline
multivariate_normal stats2random(
		const multivariate_normal& dist )
{ return dist; }

/// Uniform distribution mapping for random number generation
inline
std::uniform_real_distribution<double> stats2random(
		const std::uniform_real_distribution<double>& dist )
{ return dist; }

/// Uniform discrete distribution mapping for random number generation
inline
std::uniform_int_distribution<int> stats2random(
		const std::uniform_int_distribution<int>& dist )
{ return dist; }

/// Bernoulli distribution mapping for random number generation
inline
std::bernoulli_distribution stats2random(
		const std::bernoulli_distribution& dist )
{ return dist; }

}

/// Initializes the engine with a random seed (obtained from std::time()).
void randomize( unsigned int seed = std::time(0) );

/// Obtain a random sample from any specified distribution.
template<class T>
auto rand( T distribution ) -> decltype(_stats::stats2random(distribution)(_stats::engine))
{
	return _stats::stats2random ( distribution )( _stats::engine );
}

/// Obtain a vector of random samples from any specified distribution.
template<class T>
auto rand( T distribution, int M )
	-> vec<decltype(_stats::stats2random(distribution)(_stats::engine))>
{
	typedef decltype(_stats::stats2random(distribution)(_stats::engine)) R;
	vec<R> rv(M);
	for( int m = 0; m < M; m++ )
		rv(m) = _stats::stats2random(distribution)(_stats::engine);
	return rv;
}

/// Obtain a matrix of random samples from any specified distribution.
template<class T>
auto rand( T distribution, int M, int N )
	-> mat<decltype(_stats::stats2random(distribution)(_stats::engine))>
{
	typedef decltype(_stats::stats2random(distribution)(_stats::engine)) R;
	mat<R> rv(M,N);
	for( int m = 0; m < M; m++ )
		for( int n = 0; n < N; n++ )
			rv(m,n) = _stats::stats2random(distribution)(_stats::engine);
	return rv;
}

/// Obtain a random sample from the normal(0,1) distribution.
inline double randn() { return rand( normal(0,1) ); }

/// Obtain a random vector from the normal(0,1) distribution.
inline rvec randn( int M ) { return rand( normal(0,1), M ); }

/// Obtain a random matrix from the normal(0,1) distribution.
inline rmat randn( int M, int N ) { return rand( normal(0,1), M, N ); }

/// Obtain a random sample from the uniform(0,1) distribution.
inline double randu() { return rand( uniform(0,1) ); }

/// Obtain a random vector from the uniform(0,1) distribution.
inline rvec randu( int M ) { return rand( uniform(0,1), M ); }

/// Obtain a random matrix from the uniform(0,1) distribution.
inline rmat randu( int M, int N ) { return rand( uniform(0,1), M, N ); }
/// @}

//------------------------------------------------------------------------------
/// @name Analysis of statistical distributions
using boost::math::pdf;
using boost::math::cdf;
using boost::math::complement;
using boost::math::mean;
using boost::math::standard_deviation;
using boost::math::variance;
/// @}

/// @}
}
#endif
