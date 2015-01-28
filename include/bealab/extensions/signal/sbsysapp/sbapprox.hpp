// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/extensions/signal/sbsysapp/sbapprox.hpp
/// Base classes for subband approximation problems

#ifndef _BEALAB_EXTENSIONS_SIGNAL_SBSYSAPP_SBAPPROX_
#define	_BEALAB_EXTENSIONS_SIGNAL_SBSYSAPP_SBAPPROX_

#include <bealab/extensions/signal/sbsysapp/sbconfig.hpp>
#include <bealab/extensions/signal/sparse.hpp>

namespace bealab
{
namespace signal
{
namespace sbsysapp
{
/// @defgroup sbsysapp_sbapprox Subband approximation
/// Base classes for subband approximation problems
/// @{

/// Base for all subband approximation problems.
/// It has a subband configuration and a pool of subband indexes.
/// The subband model is formed with indexes from the pool.
class sbapprox : public sbconfig, public sbindex_pool {

protected:

	mat<cseq> G;																///< polyphase matrix of the target

public:

	using sbconfig::M;															// Use the number of subbands from the subband configuration

	/// Constructor
	sbapprox( int M_, int D_ ) : sbconfig(M_,D_), sbindex_pool(M_) {}

	/// Virtual destructor
	virtual
	~sbapprox() {}

	/// Set target (polyphase)
	virtual
	void set_target( const mat<cseq>& G_ )
	{
		assert( (int)G_.size1() == D && (int)G_.size2() == D );
		G = G_;
	}

	/// Set target (impulse response)
	void set_target( const vec<cseq>& g )
	{
		set_target( fb2pp_system(g) );
	}

	/// Set target (single complex impulse response)
	void set_target( const cseq& g0 )
	{
		vec<cseq> g(D);
		for( int d = 0; d < D; d++ )
			g(d) = g0;
		set_target( g );
	}

	/// Set target (single real impulse response)
	void set_target( const rseq& g0 )
	{
		set_target( cseq(g0) );
	}
};

/// Base for all subband approximation problems in the frequency domain.
class fdcriterion : public virtual sbapprox {

	mat<sequence<vec<cvec>>> real_atoms_sbm;											///< Table with pre-computed real atoms in the frequency domain
	mat<sequence<vec<cvec>>> imag_atoms_sbm;											///< Table with pre-computed imaginary atoms in the frequency domain
	mat<cseq> FASWA;															///< Pre-computed matrix
	mat<cseq> WSH;																///< Pre-computed matrix

	/// Get a subband model atom
	vec<cvec> get_atom_sbm_( const sbindex& sbi );

protected:

	const int Ndtft;															///< Number of points used for DTFT
	rvec wf2;																	///< Spectral weight
	rseq w;																		///< Impulse response of the square-root of the spectral weight

	/// Reset the table of subband atoms
	void reset_sbatoms();

	/// Precomputing for filterbank atoms
	void precompute_fbatoms();

	/// Get a subband model atom, and fill the table if necessary
	vec<cvec>& get_atom_sbm( const sbindex& sbi );

	/// Get an analysis window atom
	vec<cvec> get_atom_analysis( int idx );

	/// Get a synthesis window atom
	vec<cvec> get_atom_synthesis( int idx );

public:

	using sbapprox::set_target;

	vec<cvec> gf;																///< Frequency representation of the target
	function<cvec(const cseq&)> dtft;											///< N-point DTFT functor

	/// Constructor
	fdcriterion( int M, int D, int N );

	/// Virtual function override
	void set_target( const mat<cseq>& G ) override;

	/// Set the analysis filterbank
	void set_analysis_fb( const cseq& h0 );

	/// Set the synthesis filterbank
	void set_synthesis_fb( const cseq& f0 );

	/// Set the spectral weight
	void set_spectral_weight( const rseq& w_ );

	/// Objective function
	virtual double errfun()=0;
};

/// Base for frequency-domain problems using logarithmic amplitude
class logarithmic : public fdcriterion {

	/// Gradient structure
	struct gradient_t {
		rvec h0;																///< Gradient w.r.t. the analysis filterbank
		rvec f0;																///< Gradient w.r.t. the synthesis filterbank
		rvec sbm;																///< Gradient w.r.t. the subband model
	};

	figure fig_amp;																///< Figure for amplitude trace plotting
	figure fig_phase;															///< Figure for phase trace plotting

	/// Auxiliary function to compute the gradient given the atoms.
	virtual
	rvec gradient_aux( const vec<vec<cvec>>& atoms );

public:

	function<cvec(const cvec&)> logspec;										///< Log-amplitude spectrum functor

	/// Constructor
	logarithmic( int M, int D, int N );

	/// Virtual function override
	double errfun() override;

	/// Gradient (subband model)
	rvec gradient_sbm();

	/// Gradient (complete)
	gradient_t gradient();

	/// Plot function
	void plotfun();
};

/// Base for frequency-domain problems using logarithmic amplitude
/// and ignoring the phase
class logarithmic_nophase : public logarithmic {

	/// Virtual function override
	rvec gradient_aux( const vec<vec<cvec>>& atoms ) override;

public:

	/// Constructor
	logarithmic_nophase( int M, int D, int N );
};

/// @}
}
}
}
#endif
