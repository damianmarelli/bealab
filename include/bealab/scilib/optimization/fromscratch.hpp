/// @file bealab/scilib/optimization/fromscratch.hpp
/// Optimization routines developed from scratch.

#ifndef _BEALAB_OPTIMIZATION_FROMSCRATCH_
#define	_BEALAB_OPTIMIZATION_FROMSCRATCH_

#include <bealab/scilib/optimization/base.hpp>

namespace bealab
{
namespace optimization
{
/// @defgroup optimization-fromscratch From scratch
/// Optimization routines developed from scratch.
/// @{

/// Quasi-Newton method
class quasinewton : public base {

	double f0;																	///< Value of the objective function after the previous line search.
	double f;																	///< Value of the objective function after the current line search.

	/// BFGS approximation to the inverse Hessian
	rmat inverse_hessian( const rmat& H, const rvec& x, const rvec& x0,
			const rvec& g, const rvec& g0 );

	/// Line search algorithm
	rvec linesearch( const rvec& x0, const rvec& g, const rmat& B );

	/// Stopping condition for the centering algorithm
	bool stopping_condition( const rvec& x, const rvec& x0, const rvec& g );

public:

	using base::base;

	double alpha = 0;															///< Parameter alpha from [Boyd & Vandenberghe, 2004, p. 464]
	double beta  = 0.5;															///< Parameter beta from [Boyd & Vandenberghe, 2004, p. 464]

	/// Optimization routine
	rvec optimize() override;
};

/// Barrier interior-point method
class barrier : public base {

	base& cen;																	///< Centering algorithm
	double t;																	///< Centering parameter

	/// Barrier function
	double barrier_function( const rvec& x );

	/// Gradient of the barrier function
	rvec barrier_gradient( const rvec& x );

	/// Unconstrained optimization for a given path parameter 't'.
	rvec centering( const rvec& wstart );

public:

	double stop_duality_gap = 0;												///< Stopping condition

	/// Constructor
	barrier( int dim, base& c ) : base(dim), cen(c){}

	/// Optimization routine
	rvec optimize() override;
};

/// @}
}
}
#endif
