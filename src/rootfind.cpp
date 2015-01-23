// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

#include <bealab/core/gsl.hpp>
#include <bealab/scilib/optimization.hpp>
#include <bealab/scilib/calculus.hpp>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multimin.h>

namespace bealab
{

double fzero( function<double(double)> fun, double x_lo, double x_hi,
		int max_iter, double epsabs, double epsrel )
{
	// Allocate
	gsl_function F = { gsl::sfunction_proxy, &fun };
	const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
	gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set( s, &F, x_lo, x_hi );

	// Main loop
	int status;
	double rv;
	int iter  = 0;
	do {
		iter++;
		status = gsl_root_fsolver_iterate(s);
		rv     = gsl_root_fsolver_root(s);
		x_lo   = gsl_root_fsolver_x_lower(s);
		x_hi   = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval( x_lo, x_hi, epsabs, epsrel );
	} while (status == GSL_CONTINUE && iter < max_iter);

	if( iter == max_iter )
		warning("bealab::fzero - maximum number of iterations reached");

	// Free
	gsl_root_fsolver_free (s);

	return rv;
}

}
