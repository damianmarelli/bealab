/// @file bealab/scilib/calculus.hpp
/// Numeric integration and differentiation.

#ifndef _BEALAB_CALCULUS_
#define	_BEALAB_CALCULUS_

#include <bealab/core/blas.hpp>
#include <bealab/scilib/functors.hpp>

namespace bealab
{
/// @defgroup calculus Calculus
/// Numeric integration and differentiation.
/// @{

/// @name Integration

/// Integration parameters
struct integration_parameters {

	enum
	{
		dflt    = 0,
		gauss15 = 1,      // 15 point Gauss-Kronrod rule
		gauss21 = 2,      // 21 point Gauss-Kronrod rule
		gauss31 = 3,      // 31 point Gauss-Kronrod rule
		gauss41 = 4,      // 41 point Gauss-Kronrod rule
		gauss51 = 5,      // 51 point Gauss-Kronrod rule
		gauss61 = 6,      // 61 point Gauss-Kronrod rule
		montecarlo_plain = 7,
		montecarlo_miser = 8,
		montecarlo_vegas = 9
	} key;

	int subintervals;
	double epsabs;
	double epsrel;
	size_t ncalls;
	integration_parameters()
	{
		key          = dflt;
		subintervals = 10000;
		epsabs       = 1e-6;
		epsrel       = 1e-6;
		ncalls       = 1e3;
	}
};

/// Exception thrown by integral functions
class integral_error : public exception {};

/// Numeric integration of a function R->R.
double integral( const function<double(double)>&, double, double,
		integration_parameters=integration_parameters() );

/// Numeric integration of a function R->C.
complex integral( const function<complex(double)>&, double, double,
		integration_parameters=integration_parameters() );

/// Numeric integration of a function R^n->R.
double integral( const function<double(const rvec&)>&, const rvec&, const rvec&,
		integration_parameters=integration_parameters() );

/// Numeric integration of a function R^n->C.
complex integral( const function<complex(const rvec&)>&, const rvec&, const rvec&,
		integration_parameters=integration_parameters() );
/// @}

/// @name Derivation

/// First derivative of a functor at a point
double derivative( const function<double(double)>&, double );

/// N-th derivative of a functor at a point
double derivative( const function<double(double)>&, int, double );

/// N-th derivative of a functor.
auto derivative( const function<double(double)>&, int=1 ) -> function<double(double)>;

/// Gradient of a functor at a point
rvec gradient( const function<double(const rvec&)>&, const rvec& );

/// Gradient of a functor
auto gradient( const function<double(const rvec&)>& ) -> function<rvec(const rvec&)>;

/// Hessian of a functor at a point
rmat hessian( const function<double(const rvec&)>&, const rvec& );

/// Hessian of a functor
auto hessian( const function<double(const rvec&)>& ) -> function<rmat(const rvec&)>;

/// Jacobian of a functor at a point
rmat jacobian( const function<rvec(const rvec&)>&, const rvec& );

/// Jacobian of a functor
auto jacobian( const function<rvec(const rvec&)>& ) -> function<rmat(const rvec&)>;

/// @}
/// @}
}
#endif
