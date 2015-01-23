// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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
