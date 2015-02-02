/// @file bealab/scilib/optimization/ipopt.hpp
/// Wrapping of the primal-dual optimization routine from the Ipopt library.

#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOIPOPT

#ifndef _BEALAB_OPTIMIZATION_IPOPT_
#define	_BEALAB_OPTIMIZATION_IPOPT_

#include <bealab/scilib/optimization/base.hpp>

namespace bealab
{
namespace optimization
{
/// @defgroup optimization-ipopt Ipopt-based
/// Wrapping of the primal-dual optimization routine from the Ipopt library.
/// @{

/// Primal-dual method
class primal_dual : public base {

	/// Moves the data pointed by px to an rvec
	static
	rvec ptr2vec( int N, double* px );

	/// Moves the data into an rvec into the memory pointed by px
	static
	void vec2ptr( const rvec& x, double* px );

	/// Moves the data into an rmat into the memory pointed by px
	static
	void mat2ptr( const rmat& x, double* px );

	/// Objective function proxy
	static
	int eval_f( int n, double* px, int new_x, double* obj_value,
			void* user_data );

	/// Proxy for the function to compute the gradient of the objective function
	static
	int eval_grad_f( int n, double* px, int new_x, double* grad_f,
			void* user_data );

	/// Constraint function proxy
	static
	int eval_g( int n, double* px, int new_x, int m, double* g,
			void* user_data );

	/// Proxy for the function to compute the Jacobian of the objective function
	static
	int eval_jac_g( int J, double* px, int new_x, int I, int nele_jac,
			int* iRow, int *jCol, double* values, void* user_data );

	/// Proxy for computing the Hessian of the Lagrangian function
	/// For the moment it does nothing
	static
	int eval_h( int n, double* x, int new_x, double obj_factor, int m,
			double* lambda, int new_lambda, int nele_hess, int* iRow,int* jCol,
			double* values, void* user_data )
	{
		return false;
	}

public:

	/// Constructor
	primal_dual( int dim ) : base(dim) {}

	/// Optimization rutine
	rvec optimize() override;
};

/// @}
}
}
#endif
#endif
