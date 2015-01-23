// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/scilib/optimization/nlopt.hpp
/// Wrapping of optimization routines from the NLopt library.

#ifndef _BEALAB_OPTIMIZATION_NLOPT_
#define	_BEALAB_OPTIMIZATION_NLOPT_

#include <bealab/scilib/optimization/base.hpp>

namespace bealab
{
namespace optimization
{
/// @defgroup optimization-nlopt NLopt-based
/// Wrapping of optimization routines from the NLopt library.
/// @{

/// Base class for all NLopt-based algorithms
class nlopt : public base {

	/// NLopt proxy for scalar functions
	static
	double fungrad_proxy( unsigned N, const double* px, double* pgrad,
			void* pfunctors );

	/// NLopt proxy for vector functions
	static
	void vfungrad_proxy( unsigned M, double* presult, unsigned N,
			const double* px, double* pgrad, void* pfunctors );

protected:

	void* pproblem;																///< Pointer to the NLopt problem

	/// Constructor
	nlopt( int dim );

	/// Destructor
	~nlopt();

public:

	/// Optimize
	rvec optimize() override;

	void* get_pproblem() const { return pproblem; };
};

/// Nelder-Mead simplex method
class neldermead : public nlopt { public: neldermead( int dim ); };

/// Limited-memory BFGS method
class lbfgs : public nlopt { public: lbfgs( int dim ); };

/// Truncated Newton method
class tnewton : public nlopt { public: tnewton( int dim ); };

/// Truncated Newton method with preconditioning
class tnewton_precond : public nlopt { public: tnewton_precond( int dim ); };

/// Truncated Newton method with restart
class tnewton_restart : public nlopt { public: tnewton_restart( int dim ); };

/// Truncated Newton method with preconditioning and restart
class tnewton_precond_restart : public nlopt { public: tnewton_precond_restart( int dim ); };

/// Shifted limited-memory variable-metric rank-1 method
class var1 : public nlopt { public: var1( int dim ); };

/// Shifted limited-memory variable-metric rank-2 method
class var2 : public nlopt { public: var2( int dim ); };

/// Method of moving asymptotes
class mma : public nlopt { public: mma( int dim ); };

/// Sequential quadratic programming method
class sqp : public nlopt { public: sqp( int dim ); };

/// Augmented Lagrangian algorithm
class auglag : public nlopt { public: auglag( int dim, const nlopt& opt ); };

/// @}
}
}
#endif
