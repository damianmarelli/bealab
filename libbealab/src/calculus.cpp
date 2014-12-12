#include <gsl/gsl_errno.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <bealab/core/gsl.hpp>
#include <bealab/scilib/fourier.hpp>
#include <bealab/scilib/calculus.hpp>

namespace bealab
{
//-----------------------------------------------------------------------------------------------
// Integration
//-----------------------------------------------------------------------------------------------
double integral( const function<double(double)>& fun, double from, double to,
		integration_parameters ip )
{
	// Turn off error handling XXX Put this an initialization
	gsl_set_error_handler_off();

	// Allocate workspace
	gsl_integration_workspace * w =
			gsl_integration_workspace_alloc(ip.subintervals);

	// Call parameters
	gsl_function F = { _gsl::sfunction_proxy, &const_cast<function<double(double)>&>(fun) };
	int key;
	switch( ip.key ) {
	case ip.gauss15:
		key = GSL_INTEG_GAUSS15;
		break;
	case ip.gauss21:
	default:
		key = GSL_INTEG_GAUSS21;
		break;
	case ip.gauss31:
		key = GSL_INTEG_GAUSS31;
		break;
	case ip.gauss41:
		key = GSL_INTEG_GAUSS41;
		break;
	case ip.gauss51:
		key = GSL_INTEG_GAUSS51;
		break;
	case ip.gauss61:
		key = GSL_INTEG_GAUSS61;
		break;
	}

	// Return parameters
	double	result, error;
	int		retv;
//	size_t	nevals;

	if( from == to  )
		return 0;
	else if( from == -inf && to == inf )
		retv = gsl_integration_qagi( &F,
				ip.epsabs, ip.epsrel, ip.subintervals, w, &result, &error );
	else if( from == -inf && to != inf )
		retv = gsl_integration_qagil( &F, to,
				ip.epsabs, ip.epsrel, ip.subintervals, w, &result, &error );
	else if( from != -inf && to == inf )
		retv = gsl_integration_qagiu( &F, from,
				ip.epsabs, ip.epsrel, ip.subintervals, w, &result, &error );
	else
		retv = gsl_integration_qag( &F, from, to,
				ip.epsabs, ip.epsrel, ip.subintervals, key, w, &result, &error );
//		retv = gsl_integration_qng( &F, from, to,
//				epsabs, epsrel, &result, &error, &nevals );

	// free workspace
	gsl_integration_workspace_free(w);

	// analyze return value
	if(retv != 0)
		throw integral_error();
//	if(retv != 0) {
//		std::cerr << "Some problem in bealab::integral()" << std::endl;
//		abort();
//	}

//	// return error
//	if( perror != NULL )
//		*perror = error;

	return result;
}

complex integral( const function<complex(double)>& fun, double from, double to,
		integration_parameters ip )
{
	// Ral and imaginary functors
	function<double(double)> fun_real = [&fun]( double x )
	{
		return real( fun(x) );
	};
	function<double(double)> fun_imag = [&fun]( double x )
	{
		return imag( fun(x) );
	};

	return integral( fun_real, from, to, ip ) + i * integral( fun_imag, from, to, ip );
}

////XXX Replace it by a multidimensional integration method
//double integral( const function<double(const rvec&)>& fun, const rvec& a, const rvec& b,
//		integration_parameters ip )
//{
//	// Number of variables
//	assert( a.size() == b.size() );
//	int N = a.size();
//
//	// Partial integrals
//	typedef function<double(const rvec&)> VF;
//	VF partial_integral[N+1];
//	partial_integral[0] = fun;
//
//	// Form the partial integrals
//	for( int n = 0; n < N; n++ ) {
//
//		double an      = a(n);
//		double bn      = b(n);
//		VF& pintegraln = partial_integral[n];
//
//		partial_integral[n+1] = [&pintegraln,an,bn,&ip]( const rvec& x ) {
//
//			// Function of the first variable
//			function<double(double)> fn = [&](double y) {
//				rvec z = { {y}, x };
//				return pintegraln( z );
//			};
//
//			return integral( fn, an, bn, ip );
//		};
//	}
//
//	return partial_integral[N]( rvec() );
//}

double integral( const function<double(const rvec&)>& fun, const rvec& xl, const rvec& xu,
		integration_parameters ip )
{
	assert( xl.size() == xu.size() );

	// Parameters
	double result, err;
	size_t dim = xl.size();
	gsl_monte_function G = { _gsl::vfunction_proxy, dim,
			&const_cast<function<double(const rvec&)>&>(fun) };

	// Prepare environment
	gsl_rng* r = gsl_rng_alloc( gsl_rng_default );
	gsl_rng_env_setup();

	// Call the integration routine
	switch( ip.key ) {
	case ip.montecarlo_plain:
	default: {
		gsl_monte_plain_state *s = gsl_monte_plain_alloc( dim );
		gsl_monte_plain_integrate( &G, &xl(0), &xu(0), dim, ip.ncalls, r, s, &result, &err );
		gsl_monte_plain_free( s );
		break;
	}
	case ip.montecarlo_miser: {
		gsl_monte_miser_state *s = gsl_monte_miser_alloc( dim );
		gsl_monte_miser_integrate( &G, &xl(0), &xu(0), dim, ip.ncalls, r, s, &result, &err );
		gsl_monte_miser_free( s );
		break;
	}
	case ip.montecarlo_vegas: {
		rvec xl_ = xl;
		rvec xu_ = xu;
		gsl_monte_vegas_state *s = gsl_monte_vegas_alloc( dim );
		gsl_monte_vegas_integrate( &G, &xl_(0), &xu_(0), dim, ip.ncalls, r, s, &result, &err );
		gsl_monte_vegas_free( s );
		break;
	}
	}

	// Free environment
	gsl_rng_free( r );

	return result;
}

complex integral( const function<complex(const rvec&)>& fun, const rvec& a, const rvec& b,
		integration_parameters ip )
{
	// Ral and imaginary functors
	function<double(const rvec&)> fun_real = [&fun]( const rvec& x )
	{
		return real( fun(x) );
	};
	function<double(const rvec&)> fun_imag = [&fun]( const rvec& x )
	{
		return imag( fun(x) );
	};

	return integral( fun_real, a, b, ip ) + i * integral( fun_imag, a, b, ip );
}

//double integral( const function<double(int)>& fun, int from, int to,
//		integration_parameters ip )
//{
//	// Detect inverted ranges
//	int a    = min( from ,to );
//	int b    = max( from ,to );
//	int sign = (from<=to) ? 1 : -1;
//
//	// Root-integrand functor
//	function<double(int)> rintegrand = [&fun,a,b]( int t ) -> double
//	{
//		return (t<a || t > b) ? 0 : sqrt( fun(t) );
//	};
//
//	// Dual functor
//	function<double(double)> dualfun = [&rintegrand]( double w )
//	{
//		return pow( dtft(rintegrand)( w ), 2 );
//	};
//
//	return sign * integral( dualful, -inf, inf );
//}

//-----------------------------------------------------------------------------------------------
// Derivatives
//-----------------------------------------------------------------------------------------------
double derivative( const function<double(double)>& fun, double x )
{
	return derivative( fun, 1, x );
}

double derivative( const function<double(double)>& fun, int n, double x )
{
	// Prepare the function to derive
	function<double(double)> fund;
	if( n == 1 )
		fund = fun;
	else
		fund = derivative( fun, n-1 );

	// Call parameters
	gsl_function F = { _gsl::sfunction_proxy, &fund };
	double result, abserr;
	double stepsize = 1e-4*abs(x)+1e-4;

	// Derive
	gsl_deriv_central (&F, x, stepsize, &result, &abserr);

	return result;
}

auto derivative( const function<double(double)>& f, int n ) -> function<double(double)>
{
	return [f,n](double x) { return derivative( f, n, x ); };
}

rvec gradient( const function<double(const rvec&)>& f, const rvec& x )
{
	int N = x.size();
	rvec g(N);
	for( int n = 0; n < N; n++ ) {
		auto fn = [&](double xn) {
			rvec x_ = x;
			x_(n) = xn;
			return f(x_); };
		g(n) = derivative( fn, x(n) );
	}

	return g;
}

auto gradient( const function<double(const rvec&)>& f ) -> function<rvec(const rvec&)>
{
	return [f](const rvec &x) { return gradient( f, x ); };
}

rmat hessian( const function<double(const rvec&)>& f, const rvec &x )
{
	auto g = [&](const rvec &x_) { return gradient(f,x_); };
	rmat H = jacobian( g, x );
	return (H + trans(H)) / 2;
}

auto hessian( const function<double(const rvec&)>& f ) -> function<rmat(const rvec&)>
{
	return [f](const rvec &x) { return hessian( f, x ); };
}

rmat jacobian( const function<rvec(const rvec&)>& f, const rvec &x )
{
	rvec y = f(x);
	int M  = y.size();
	int N  = x.size();
	rmat J(M,N);
	for( int m = 0; m < M; m++ ) {
		auto fm = [&](const rvec &x_) { return f(x_)(m); };
		J.row( m ) = gradient( fm, x );
	}

	return J;
}

auto jacobian( const function<rvec(const rvec&)>& f ) -> function<rmat(const rvec&)>
{
	return [f](const rvec &x) { return jacobian( f, x ); };
}

}
