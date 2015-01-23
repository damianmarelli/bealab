// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/extensions/signal/sbsysapp/sbrelax.hpp
/// Subband system approximation using relaxation methods.

#ifndef _BEALAB_SBSYSAPP_SBRELAX_
#define	_BEALAB_SBSYSAPP_SBRELAX_

#include <bealab/extensions/signal/sbsysapp/sbgreedy.hpp>

namespace bealab
{
namespace signal
{
namespace sbsysapp
{
/// @defgroup sbsysapp_sbrelax Relaxation methods
/// Subband system approximation using relaxation methods.
/// @{

/// Subband approximation using Hyder-Mahata's relaxation.
template<class criterion>
class relaxation :
	public criterion,
	public sparse::shrinking_gaussian {

	timer tim;

	figure fig = figure("Coefficient distribution");							///< Figure for plotting the coefficient's distribution

protected:

		/// Plotting function
		void plotfun()
		{
			if( tim.elapsed() < 1 )
				return;

			// Reset timer
			tim.reset();

			// Plot the base
//			criterion::plotfun();

			// Local plot
			rvec coeffs = abs(this->sbm.sbvalues);
			sort( coeffs.begin(), coeffs.end() );
			coeffs = flip(coeffs);
//			coeffs = coeffs / coeffs(0);
			fig.clear().logy().plot( coeffs );
		}

//		/// Obtain a guess using the nlpursuit method
//		sbmodel get_guess()
//		{
//			// Run a nlpursuit subproblem
//			nlpursuit<criterion> sbapp( this->M, this->D, this->Ndtft );
//			sbapp.set_analysis_fb( this->h0() );
//			sbapp.set_synthesis_fb( this->f0() );
//			sbapp.set_sbindexes( this->sbindexes );
//			sbapp.set_target( this->G );
//			sbapp.set_spectral_weight( this->w );
//			sbapp.tolerance = this->tolerance;
//			sbapp.trace     = trace;
//			sbapp.plot      = plot;
//			sbapp.approximate();
//
//			return sbapp.sbm;
//		}

		/// Tune the result
		void tune()
		{
			optimization::var2 opt( this->sbm.sbvalues.size() );
			opt.set_objective(
					[this](const rvec& x)
					{
						this->sbm.sbvalues = x;
						return dynamic_cast<criterion*>(this)->errfun();
					},
					[this](const rvec& x)
					{
						this->sbm.sbvalues = x;
						return this->gradient_sbm();
					} );
			opt.guess          = this->sbm.sbvalues;
			opt.trace          = this->trace;
			this->sbm.sbvalues = opt.optimize();
		}

public:

	int max_sparsity = 0;														///< Dummy member to present a unified interface with sbsysapp::ompursuit and sbsysapp::nlpursuit

	/// Constructor
	relaxation( int M, int D, int N ) :
		sbapprox(M,D), criterion(M,D,N)
	{
		// Approximation error function
		errfun = [this]( const rvec& x )
		{
			this->sbm.sbindexes = this->sbindexes;
			this->sbm.sbvalues  = x;
			double err = this->criterion::errfun();
			if( this->trace >= 2 )
				cout << "subband::relaxation - Error = " << err << endl;
			if(this->plot)
				this->plotfun();
			return err;
		};

		// Gradient of the approximation error function
		errfun_gradient = [this]( const rvec& x )
		{
			this->sbm.sbindexes = this->sbindexes;
			this->sbm.sbvalues  = x;
			return this->gradient().sbm;
		};
	}

	/// Set the initialization
	void set_guess( const sbmodel& gss )
	{
		this->sbindexes = gss.sbindexes;
		guess           = gss.sbvalues;
	}

	/// Virtual function override
	void approximate() override
	{
//		// Initialize the support of the subband indexes and obtain a guess
//		// for the coefficients
//		sbmodel gss     = get_guess();
//		this->sbindexes = gss.sbindexes;
//		guess           = gss.sbvalues;

		// Sparsification using the Hyder-Mahata  method
		shrinking_gaussian::approximate();

		// Form the subband model with the sparsified parameters
		this->sbm.sbindexes = this->sbindexes(indirect(indexes));
		this->sbm.sbvalues  = coefficients;

		// Tune the result
		tune();
	}
};

/// @}
}
}
}
#endif
