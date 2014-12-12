#include <bealab/core/blas.hpp>

namespace bealab
{

rvec vslice( double from, double skip, double N )
{
//	int N = abs( floor( (to-from) / skip ) + 1 );
	rvec x(N);
	for( int n = 0; n < N; n++ )
		x(n) = from + n*skip;
	return x;
}

ivec vrange( int from, int to_plus_1 )
{
	return vslice( from, 1, to_plus_1-from );
}

rvec linspace( double from, double to, int N )
{
//	if( N < 0 ) {
//		assert(  mod(to-from,1) == 0 );
//		N = abs(to-from) + 1;
//	}

	// Case N == 0
	if( N == 0 )
		return rvec();

	// Case from == to
	if( from == to )
		return from * ones(N);

	// Case N == 1 and from ~= to
	if( N == 1 )
		error("linspace() - At least two points need to be generated");

	// Normal case: N > 1 and from ~= to
	double skip = (to-from) / (N-1);
	rvec x(N);
	for( int n = 0; n < N; n++ )
		x(n) = from + n*skip;
	return x;
}

rvec logspace( double a, double b, int N )
{
	double alog = log(a);
	double blog = log(b);
	rvec xlog   = linspace( alog, blog, N );
	return exp(xlog);
}

}
