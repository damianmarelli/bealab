/// @file bealab/scilib/optimization/gsl.hpp
/// Wrapping of optimization routines from the GSL library.

#ifndef _BEALAB_OPTIMIZATION_GSL_
#define	_BEALAB_OPTIMIZATION_GSL_

#include <bealab/scilib/optimization/base.hpp>

namespace bealab
{
namespace optimization
{
/// @defgroup optimization-gsl GSL-based
/// Wrapping of optimization routines from the GSL library.
/// @{

/// Base class for all GSL-based algorithms
class gsl : public base {
protected:

	void *pproblem;																///< Pointer to the GSL problem

	/// Functors passed as parameters to the proxies
	struct functors {
		function<double(const rvec&)>        fun;
		function<rvec (const rvec&)>         grad;
		function<double(const rvec&, rvec&)> fungrad;
	};

	/// Proxy for the objective function
	static
	double objective_proxy( const void* x, void* par );

	/// Proxy for the gradient
	static
	void gradient_proxy( const void* x, void* par, void* gr );

	/// Proxy for the joint objective-gradient function
	static
	void objgrad_proxy( const void* x, void* par, double* f, void* gr );

	/// @name Controlling the stopping condition
	int status;
	double fx, fx0;
	rvec x, x0;
	rvec grad;
	int Nfeval;
	bool check_stopping_condition();
	/// @}

	/// name Constructor
	gsl( int dim ) : base(dim) {}
};

/// Base class for all GSL-based derivative-free algorithms
class gsl_noder : public gsl {
public:

	/// Constructor
	gsl_noder( int dim ) : gsl(dim){}

	/// Do the optimization
	rvec optimize() override;
};

/// Base class for all GSL-based derivative-based algorithms
class gsl_der : public gsl {
public:

	bool test_gradient = false;													///< If true, the gradient is tested by comparing it with the numerical one.

	/// Constructor
	gsl_der( int dim ) : gsl(dim){}

	/// Do the optimization
	rvec optimize() override;
};

/// Simplex method
class simplex : public gsl_noder { public: simplex( int dim ); };

/// Conjugate-gradient method from Fletcher-Reeves
class conjugate_fr : public gsl_der { public: conjugate_fr( int dim ); };

/// Conjugate-gradient method from Polak-Ribiere
class conjugate_pr : public gsl_der { public: conjugate_pr( int dim ); };

/// Broyden-Fletcher-Goldfarb-Shanno (BFGS) method
class bfgs : public gsl_der { public: bfgs( int dim ); };

/// @}
}
}
#endif
