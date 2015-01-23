// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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
