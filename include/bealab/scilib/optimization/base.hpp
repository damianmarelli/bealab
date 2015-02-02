/// @file bealab/scilib/optimization/base.hpp
/// Base for all optimization classes.

#ifndef _BEALAB_OPTIMIZATION_BASE_
#define	_BEALAB_OPTIMIZATION_BASE_

#include <bealab/core/blas.hpp>

namespace bealab
{
/// Parameter optimization routines.
namespace optimization
{
/// @defgroup optimization-base Base
/// Base for all optimization classes.
/// @{

/// Base class for all optimization algorithms
class base {
protected:

	/// Function / gradient pair
	struct fungrad {
		std::function<double(const rvec&)> function;
		std::function<rvec(const rvec&)> gradient;
	};

	/// Function / gradient pair (vector version)
	struct vfunjac {
		std::function<rvec(const rvec&)> function;
		std::function<rmat(const rvec&)> jacobian;
	};

	/// Function / gradient pair (vector / sparse version)
	struct vsfunjac {
		std::function<rvec(const rvec&)> function;
		std::function<sparse_matrix<double>(const rvec&)> jacobian;
	};

	/// @name Objective and constraints
	fungrad objective;
	std::vector<fungrad> inequality_constraints;
	std::vector<fungrad> equality_constraints;
	std::vector<vfunjac> inequality_vconstraints;
	std::vector<vfunjac> equality_vconstraints;
	std::vector<vsfunjac> inequality_vsconstraints;
	std::vector<vsfunjac> equality_vsconstraints;
	/// @}

	/// Test the gradient of a scalar function by comparing it with the
	/// numeric one
	void test_gradient( const fungrad& fg, const rvec& x0 );

	/// Test the Jacobian of a vector function by comparing it with the
	/// numeric one
	void test_jacobian( const vfunjac& fj, const rvec& x0 );

public:

	/// Stop codes
	enum stopcode { none, error, fvalue, fincrement, fincrement_relative,
		fincrement_absolute, xincrement, xincrement_relative,
		xincrement_absolute, gradient, feval, time, force };

	/// @name Optimization parameters
	int dim;																	///< Dimension of the parameter vector
	int trace = 0;
	rvec guess;
	rvec lower_bounds;
	rvec upper_bounds;
	/// @}

	/// @name Stopping conditions
	double stop_gradient             = 0;
	int stop_feval                   = imax;
	double stop_time                 = inf;
	double stop_fvalue               = -inf;
	double stop_fincrement_relative  = 0;
	double stop_fincrement_absolute  = 0;
	double stop_xincrement_relative  = 0;
	rvec  stop_xincrement_absolute;
	double stop_constraint_tolerance = 1e-8;
	bool stop_force                  = false;
	/// @}

	/// @name Result
	rvec solution;
	double minimum      = nan;
	stopcode stopreason = none;
	/// @}

	/// Callback function
	function<void(const rvec&)> callback_function = 0;

	/// @name Constructor & destructor
	base( int dim_ )
	{
		dim                      = dim_;
		guess                    = zeros(dim);
		lower_bounds             = -inf * ones(dim);
		upper_bounds             =  inf * ones(dim);
		stop_xincrement_absolute = zeros(dim);
	}

	virtual
	~base() {};
	/// @}

	/// @name Set functions
	void set_objective( const function<double(const rvec&)>& fun,
						const function<rvec(const rvec&)>& grad=0 );

	void add_inequality_constraint( const function<double(const rvec&)>& fun,
						 const function<rvec(const rvec&)>& grad=0 );

	void add_equality_constraint( const function<double(const rvec&)>& fun,
						 const function<rvec(const rvec&)>& grad=0 );

	void add_inequality_vconstraint( const function<rvec(const rvec&)>& fun,
						  const function<rmat(const rvec&)>& jacob=0 );

	void add_equality_vconstraint( const function<rvec(const rvec&)>& fun,
						  const function<rmat(const rvec&)>& jacob=0 );

	void add_inequality_vsconstraint( const function<rvec(const rvec&)>& fun,
						  const function<sparse_matrix<double>(const rvec&)>& jacob );

	void add_equality_vsconstraint( const function<rvec(const rvec&)>& fun,
						  const function<sparse_matrix<double>(const rvec&)>& jacob );
	/// @}

	/// @name Optimize
	virtual
	rvec optimize() = 0;

	/// Test the derivatives of the objective function and constraints
	void test_derivatives( const rvec& x0 );
	/// @}
};

/// @}
}
}
#endif
