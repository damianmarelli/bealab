// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/extensions/signal/sbsysapp/sbgreedy.hpp
/// Greedy subband methods for system approximation.

#ifndef _BEALAB_EXTENSIONS_SIGNAL_SBSYSAPP_SBGREEDY_
#define	_BEALAB_EXTENSIONS_SIGNAL_SBSYSAPP_SBGREEDY_

#include <bealab/scilib/cepstrum.hpp>
#include <bealab/scilib/stats.hpp>
#include <bealab/extensions/signal/sbsysapp/sbapprox.hpp>

namespace bealab
{
namespace signal
{
namespace sbsysapp
{
/// @defgroup sbsysapp_sbgreedy Greedy methods
/// Greedy subband methods for system approximation.
/// @{

/// Base for all subband approximation problems using greedy algorithms
class greedy :
	public virtual sbapprox,
	public virtual sparse::greedy<double> {

	figure fig;																	///< Figure for trace plotting.

	/// Count the number of non-zero coefficients in a mat<cseq>
	int number_of_coefficients( const mat<cseq>& S );

protected:

	/// Virtual function override
	bool forcestop() override;

	/// Virtual function override
	void plotfun() override;

public:

	/// Constructor
	greedy( int M, int D );

	/// Virtual function override
	void approximate() override;
};

/// Kronecker subband approximation problems.
namespace kronecker
{

/// Base for Kronecker approximation problems
class linear : public virtual sbapprox {

	figure fig;																	///< Figure for trace plotting.

public:

	/// Constructor
	linear( int M, int D ) :
		sbapprox(M,D),
		fig("spectra") {}

	/// Plot function
	void plotfun( const ivec& idxs, const rvec& coeffs );
};

/// Time-domain version.
class time_domain :
	public greedy,
	public linear,
	public sparse::orthogonal_matching_pursuit<double,mat<cseq>> {

private:

	mat<cseq> HH;																///< Precomputed H*H.star()
	mat<cseq> FF;																///< Precomputed F*F.star()

	/// Callback function
	mat<cseq> synthesis_( const ivec& idxs, const rvec& coeffs );

	/// Callback function
	rvec analysis_( const mat<cseq>& G );

	/// Callback function
	double inner_product_atoms_( int i, int j );

protected:

	/// Virtual function implementation
	ivec next_indexes() override;

	/// Virtual function implementation
	void plotfun() override;

public:

	// Imports from subband
	using sbapprox::set_target;
	using sbapprox::set_analysis_fb;
	using sbapprox::set_synthesis_fb;

	/// Constructor
	time_domain( int M, int D );

	/// Virtual function override
	void set_target( const mat<cseq>& G ) override;

	/// Virtual function override
	void set_analysis_fb( const mat<cseq>& H ) override;

	/// Virtual function override
	void set_synthesis_fb( const mat<cseq>& F ) override;

	/// Virtual function override
	void set_sbindexes( const vec<sbindex>& idxs ) override;

	/// Virtual function override
	void approximate() override;
};

/// Frequency domain.
class frequency_domain :
	public greedy,
	public linear,
	public sparse::orthogonal_matching_pursuit<double,vec<cmat>> {

	mat<cvec> HHf;																///< Precomputed H*H.star()
	mat<cvec> FFf;																///< Precomputed F*F.star()
	vec<cmat> Hf;																///< Frequency version of H
	vec<cmat> Ff;																///< Frequency version of F

	/// Callback function
	vec<cmat> synthesis_( const ivec& idxs, const rvec& coeffs );

	/// Callback function
	rvec analysis_( const vec<cmat>& Gf );

	/// Callback function
	double inner_product_atoms_( int i, int j );

protected:

	/// Virtual function override
	ivec next_indexes() override;

	/// Virtual function override
	void plotfun() override;

	function<cvec(const cseq&)> dtft;
	function<cseq(const cvec&)> idtft;

public:

	function<vec<cmat>(const mat<cseq>&)> td2fd;
	function<mat<cseq>(const vec<cmat>&)> fd2td;

	// Imports from subband
	using sbapprox::set_target;
	using sbapprox::set_analysis_fb;
	using sbapprox::set_synthesis_fb;
	using sbapprox::set_sbindexes;

	/// Constructor
	frequency_domain( int M, int D, int N );

	/// Set the target
	void set_target( const vec<cmat>& Gf );

	/// Virtual function override
	void set_target( const mat<cseq>& G ) override;

	/// Set the analysis filterbank
	void set_analysis_fb( const vec<cmat>& Hf_ );

	/// Virtual function override
	void set_analysis_fb( const mat<cseq>& H ) override;

	/// Set the synthesis filterbank
	void set_synthesis_fb( const vec<cmat>& Ff_ );

	/// Virtual function override
	void set_synthesis_fb( const mat<cseq>& F ) override;

	/// Virtual function override
	void set_sbindexes( const vec<sbindex>& idxs ) override;

	/// Virtual function override
	void approximate() override;
};

}

/// Subband system approximation using orthogonal matching pursuit
template<class criterion>
class ompursuit : public kronecker::time_domain {

	criterion crit;
	rseq w;																		///< Spectral weight
	mat<rseq> W;																///< Polyphase representation of the spectral weight factor

protected:

	/// Virtual function override
	double errfun( const ivec& idxs, const rvec& coeffs ) override
	{
		// Run errfun() from kronecker::time_domain
		kronecker::time_domain::errfun( idxs, coeffs );

		// Compute the error according to criterion
		crit.sbm.sbindexes = sbindexes(indirect(idxs));
		crit.sbm.sbvalues  = coeffs;
		return crit.errfun();
	}

	/// Virtual function override
	void plotfun() override
	{
		crit.sbm.sbindexes = sbindexes(indirect(indexes));
		crit.sbm.sbvalues  = coefficients;
		crit.plotfun();
		greedy::plotfun();
	}

public:

//	using sbapprox::set_target;

	/// Constructor
	ompursuit( int M, int D, int N ) :
		sbapprox(M,D),
		kronecker::time_domain(M,D),
		crit(M,D,N)
	{
		set_spectral_weight( {1} );
	}

//	/// Virtual function override
//	void set_target( const mat<cseq>& G ) override
//	{
//		kronecker_td::set_target(G);
//		criterion::set_target(G);
//	}
//
//	/// Member function override
//	void set_analysis_fb( const rseq& h0 )
//	{
//		kronecker_td::set_analysis_fb(h0);
//		criterion::set_analysis_fb(h0);
//	}
//
//	/// Member function override
//	void set_synthesis_fb( const rseq& f0 )
//	{
//		kronecker_td::set_synthesis_fb(f0);
//		criterion::set_synthesis_fb(f0);
//	}

	/// Set the spectral weight
	void set_spectral_weight( const rseq& w_ )
	{
//		criterion::set_spectral_weight(w);
		w = w_;
		W = fb2pp_system( w, D );
	}

	/// Virtual function override
	void approximate() override
	{
		// Set the error criterion
		crit.set_analysis_fb( h0() );
		crit.set_synthesis_fb( f0() );
		crit.set_target( G );
		crit.set_spectral_weight( w );

		// Run the approximation
		mat<cseq> G_ = G;
		mat<cseq> H_ = H;
		set_target( G*W );
		set_analysis_fb( H*W );
		kronecker::time_domain::approximate();
		set_target( G_ );
		set_analysis_fb( H_ );
	}
};

/// Subband system approximation using the old orthogonal matching pursuit method
template<class criterion>
class ompursuit_old : public kronecker::frequency_domain {

	criterion crit;
	rseq w;																		///< Spectral weight
	mat<rseq> W;																///< Polyphase representation of the spectral weight factor
	vec<cmat> Rf;																///< Frequency-polyphase representation of the modified spectral weights

	/// Virtual function override
	double errfun( const ivec& idxs, const rvec& coeffs ) override
	{
		// Run errfun() from kronecker::frequency_domain
		kronecker::frequency_domain::errfun( idxs, coeffs );

		// Compute the error according to criterion
		crit.sbm.sbindexes = sbindexes(indirect(idxs));
		crit.sbm.sbvalues  = coeffs;
		return crit.errfun();
	}

	/// Virtual function override
	void plotfun() override
	{
		crit.sbm.sbindexes = sbindexes(indirect(indexes));
		crit.sbm.sbvalues  = coefficients;
		crit.plotfun();
		greedy::plotfun();
	}

	/// Update the spectral weight representing the logarithmic amplitude scale
	void update_frweight()
	{
		// Compute linear difference
		mat<cseq> Gh  = polyphase();
		mat<cseq> Gt  = G - Gh;
		vec<cmat> Gtf = td2fd( Gt );

		// Compute log difference
		vec<cvec> gf  = entrywise(dtft)( pp2fb_system(G) );
		vec<cvec> ghf = entrywise(dtft)( pp2fb_system( Gh ) );
		vec<cvec> ctf(D);
		for( int d = 0; d < D; d++ )
			ctf(d) = crit.logspec( element_div( crit.gf(d), ghf(d) ) );
		vec<cseq> ct   = entrywise( [this](const cvec& x){ return idtft(x); } )( ctf );
		mat<cseq> Ct   = fb2pp_system( ct );
		vec<cmat> Ctf  = td2fd( Ct );
		vec<cmat> Wf   = td2fd( W );
		vec<cmat> CtWf = element_prod( Ctf, Wf );

		// Compute Rf (the spectral compensation for logarithmic amplitude)
		Rf = element_prod( element_inv(Gtf), CtWf );
	}

protected:

	/// Virtual function override
	ivec next_indexes() override
	{
		ivec idxs   = orthogonal_matching_pursuit::next_indexes();
		sbindex sbi = sbindexes(idxs(0));
		if( !sbi.is_selfconjugate(M) ) {
			sbindex sbicc = sbi;
			sbicc.real_f  = !sbicc.real_f;
			idxs          = { idxs, {sbi2int(sbicc)} };
		}
		if(this->trace)
			cout << "Chosen indexes: " << sbindexes(indirect(idxs)) << endl;
		return idxs;
	}

public:

	/// Constructor
	ompursuit_old( int M, int D, int N ) :
		sbapprox(M,D),
		kronecker::frequency_domain(M,D,N),
		crit(M,D,N)
	{
		set_spectral_weight( {1} );
	}

	/// Set the spectral weight
	void set_spectral_weight( const rseq& w_ )
	{
		w = w_;
		W = fb2pp_system( w, D );
	}

	/// Virtual function override
	void approximate() override
	{
		// Set the error criterion
		crit.set_analysis_fb( h0() );
		crit.set_synthesis_fb( f0() );
		crit.set_target( G );
		crit.set_spectral_weight( w );

		// Remember the actual target and analysis FB
		mat<cseq> G_ = G;
		mat<cseq> H_ = H;

		// Main loop
		for( int i = 0;; i++ ) {

			// Remember the current subband model
			sbmodel sbm0 = sbm;

			// Approximate
			kronecker::frequency_domain::approximate();

			// Stopping condition
			if( sbm.sparsity() >= sbm0.sparsity() && i > 0 ) {
				sbm = sbm0;
				break;
			}

			// Update the problem
			update_frweight();
			vec<cmat> Hf  = td2fd( H );
			vec<cmat> Gf  = td2fd( G );
			vec<cmat> Gzf = element_prod( Gf, Rf );
			vec<cmat> Hzf = element_prod( Hf, Rf );
			set_analysis_fb( Hzf );
			set_target( Gzf );
		}

		// Reset the actual target and analysis FB
		set_analysis_fb( G_ );
		set_target( H_ );
	}
};

/// Subband system approximation using non-linear pursuit
template<class criterion>
class nlpursuit :
	public greedy,
	public criterion,
	public sparse::gradient_pursuit<mat<cseq>> {

	/// Subproblem class. Gives access to protected members from bomp_ppfd
	class subproblem_t :
		public kronecker::frequency_domain {
	public:
		subproblem_t(int M, int D, int N) :
			sbapprox(M,D), kronecker::frequency_domain(M,D,N) {}
		using kronecker::frequency_domain::optimal_coefficients;
		using kronecker::frequency_domain::synthesis;
		using kronecker::frequency_domain::next_indexes;
	};

	subproblem_t subproblem;													///< Subproblem used for obtaining the next index and the guess
	vec<cmat> Gf;																///< Frequency-polyphase representation of the target
	vec<cmat> Hf;																///< Frequency-polyphase representation of the analysis filterbank
	vec<cmat> Wf;																///< Frequency-polyphase representation of the nominal spectral weight
	vec<cmat> Rf;																///< Frequency-polyphase representation of the modified spectral weights
	bool gradientpursuit_f;														///< Flag to indicate if the next indexes are chosen using gradient pursuit

	/// Update the auxiliar variables Zf and CtWf
	void update_auxiliar_variables()
	{
//		// Compute linear difference
//		vec<cseq> gh  = pp2fb_system( approximation );
//		vec<cvec> ghf = gh.apply( dtft );
//		vec<cvec> gtf = gf - ghf;
//
//		// Compute log difference
//		vec<cvec> ctf(D);
//		for( int d = 0; d < D; d++ )
//			ctf(d) = logspec( element_div( gf(d), ghf(d) ) );
//
//		// Spectral compensation weights
//		vec<cvec> af(D);
//		for( int d = 0; d < D; d++ )
//			af(d) = element_div( ctf(d), gtf(d) );
//
//		// Average weight
//		rvec af1 = sum( abs( af ) ) / D;
//		rvec zf1 = element_prod( af1, sqrt(wf) );
//		vec<rvec> zf(D);
//		for( int d = 0; d < D; d++ )
//			zf(d) = zf1;
//		cseq z1      = idtft( zf1 );
//		mat<cseq> Z  = fb2pp_system( z1, D );
//		Zf           = sp_indexes.td2fd( Z );
//
//		// Compute the spectral compensation for logarithmic amplitude
//		if( indexes.size() == 0 )
//			Zf  = Wf;

		// Compute linear difference
		mat<cseq> Gh  = polyphase();
		mat<cseq> Gt  = G - Gh;
		vec<cmat> Gtf = subproblem.td2fd( Gt );

		// Compute log difference
		vec<cseq> gh  = pp2fb_system( Gh );
		vec<cvec> ghf = entrywise(this->dtft)( gh );
		vec<cvec> ctf(D);
		for( int d = 0; d < D; d++ )
			ctf(d) = this->logspec( element_div( this->gf(d), ghf(d) ) );
		vec<cseq> ct   = entrywise( [](const cvec& x){ return idtft(x); } )( ctf );
		mat<cseq> Ct   = fb2pp_system( ct );
		vec<cmat> Ctf  = subproblem.td2fd( Ct );
		vec<cmat> CtWf = element_prod( Ctf, Wf );

		// Compute Rf (the spectral compensation for logarithmic amplitude)
//		if( indexes.size() == 0 ) {
//			Rf = Wf;
//		}
//		else {
//			Rf = element_prod( element_inv(Gtf), CtWf );
//		}
		Rf = element_prod( element_inv(Gtf), CtWf );

		// If not well conditioned, switch to the linear weight.
		double cond = abs(sum(sum(Rf)));
		if( isnan(cond) || cond == inf )
			Rf = Wf;
	}

	/// Choose the next indexes using iterative re-weighting
	ivec next_indexes_itreweighting()
	{
		// Modify the approximation subproblem
		vec<cmat> Gzf = element_prod( Gf, Rf );
		vec<cmat> Hzf = element_prod( Hf, Rf );
		subproblem.set_analysis_fb( Hzf );
		subproblem.set_synthesis_fb( F );
		subproblem.set_target( Gzf );
		subproblem.initialization();

		// Get a new index from the subproblem
		subproblem.indexes       = indexes;
		subproblem.coefficients  = subproblem.optimal_coefficients();
		subproblem.approximation = subproblem.synthesis( subproblem.indexes, subproblem.coefficients );
		int idx                  = subproblem.next_indexes()(0);

//		// Convert it to an index of the main problem
//		SBI idx = sidx;
//		if( idx.m >= M ) {
//			idx.num_f = false;
//			idx.m    -= M;
//			idx.n    -= M;
//		}

		// Make a vector of indexes and return
		return {idx};
	}

	bool diag_is_filled( const ivec& idxs, const rvec& coeffs )
	{
		sbm.sbindexes = sbindexes(indirect(idxs));
		sbm.sbvalues  = coeffs;
		mat<cseq> S   = sbm.models().S;
		for( int m = 0; m < M; m++ )
			if( S(m,m).size() == 0 )
				return false;
		return true;
	}

	/// Virtual function override
	void set_analysis_fb( const mat<cseq>& H ) override
	{
		// Set the base
		greedy::set_analysis_fb(H);

		// Frequency-polyphase representation of the analysis filterbank
		Hf = subproblem.td2fd( H );
	}

	/// Virtual function override
	void set_synthesis_fb( const mat<cseq>& F ) override
	{
		// Set the base
		greedy::set_synthesis_fb(F);

		// Compute frequency-polypase of the dual
		vec<cmat> Ff = subproblem.td2fd( F );
	}

protected:

	/// Virtual function override
//	rvec guess() override
//	{
//		rvec coeffs;
////		if( diag_is_filled( indexes.range(0,end-1), coefficients ) ) {
////		if(false) {
//		if( (int)indexes.size() > M/2 ) {
////			cout << "Guess continue" << endl;
//			subproblem.set_target( Tf );
//			subproblem.initialization();
//			subproblem.indexes = indexes( range(indexes.size()-1,indexes.size()) );
//			rvec ncoeff        = subproblem.optimal_coefficients();
//			coeffs             = { coefficients, ncoeff };
//		}
//		else {
////			cout << "Guess restart" << endl;
//			subproblem.indexes = indexes;
//			coeffs             = subproblem.optimal_coefficients();
//		}
//		return coeffs;
//	}
	rvec guess() override
	{
//		if( indexes.size() <= 1 ) {
//			subproblem.indexes = indexes;
//			return subproblem.optimal_coefficients();
//		}
//		else
//			return { coefficients, {0} };

		// Obtain a guess
		rvec gss;
		if( indexes.size() <= 1 ) {
			subproblem.indexes = indexes;
			gss = subproblem.optimal_coefficients();
		}
		else
			gss = { coefficients, {0} };

		// Make the guess well conditioned
		for( int i = 0; i < 10; i++ ) {
			if( errfun( indexes, gss ) == inf ) {
				double amp = 1e-6 * max(abs(gss));
				gss += amp * randn(gss.size());
			}
			else
				break;
		}
		if( errfun( indexes, gss ) == inf )
			warning("Bad guess");

		return gss;
	}

	/// Virtual function override
	ivec next_indexes() override
	{
		// Update the subproblem
		update_auxiliar_variables();

		// Check the condition to switch to gradient pursuit
		if( error < error_th && gradientpursuit_f == false ) {
//		if( (int)indexes.size() >= M/2 && gradientpursuit_f == false ) {
			gradientpursuit_f = true;
			cout << "*** Switching to gradient pursuit ***" << endl;
		}

		// Choose the next index
		ivec idxs;
		if( gradientpursuit_f )
			idxs = gradient_pursuit::next_indexes();
		else
			idxs = next_indexes_itreweighting();

		// Trace
		if(this->trace)
			cout << "Chosen index: " << sbindexes(idxs(idxs.size()-1)) << endl;

		return idxs;
	}

	/// Virtual function override
	double errfun( const ivec& idxs, const rvec& coeffs ) override
	{
		sbm.sbindexes = sbindexes(indirect(idxs));
		sbm.sbvalues  = coeffs;
		return criterion::errfun();
	}

	/// Virtual function override
	rvec gradient( const rvec& coeffs, gradient_type gt ) override
	{
		if( gt == reduced ) {
			sbm.sbindexes = sbindexes(indirect(indexes));
			sbm.sbvalues  = coeffs;
		}
		else {
			sbm.sbindexes                   = sbindexes;
			sbm.sbvalues                    = zeros(sbindexes.size());
			sbm.sbvalues(indirect(indexes)) = coeffs;
		}
		return criterion::gradient_sbm();
	}

	figure fig;

	/// Virtual function override
	void plotfun() override
	{
		sbm.sbindexes = sbindexes(indirect(indexes));
		sbm.sbvalues  = coefficients;
		criterion::plotfun();
		greedy::plotfun();

//		mat<cseq> R   = subproblem.fd2td( Rf );
//		vec<cseq> r   = pp2fb_system( R );
//		vec<cvec> rf  = entrywise(this->dtft)( r );
//		fig.clear().overlap().plot( abs(rf) );
//		cin.get();
	}

public:

	// Imports from the base
	using greedy::set_target;
	using greedy::set_sbindexes;
	using criterion::set_analysis_fb;
	using criterion::set_synthesis_fb;

	double error_th;															///< Error threshold to switch to gradient pursuit

	/// Constructor
	nlpursuit( int M, int D, int N ) :
		sbapprox(M,D), greedy(M,D), criterion(M,D,N),
		subproblem(M,D,N), gradientpursuit_f(false),
		error_th(0)
	{
		set_spectral_weight( {1} );
	}

	/// Virtual function override
	void set_target( const mat<cseq>& G ) override
	{
		// Set the bases
		criterion::set_target( G );

		// Set local variables
		Gf = subproblem.td2fd( G );
	}

	/// Virtual function override
	void set_sbindexes( const vec<sbindex>& idxs ) override
	{
		// Set the bases
		greedy::set_sbindexes(idxs);
		gradient_pursuit::set_size(sbindexes.size());

		// Set the subproblem
		int N = idxs.size();
		vec<sbindex> sidxs(N);
		for( int n = 0; n < N; n++ ) {
			sidxs(n) = idxs(n);
			if( !sidxs(n).num_f ) {
				sidxs(n).num_f = true;
				sidxs(n).m    += M;
				sidxs(n).n    += M;
			}
		}
		subproblem.set_sbindexes(sidxs);
	}

	/// Member function override
	void set_spectral_weight( const rseq& w )
	{
		criterion::set_spectral_weight( w );
		mat<cseq> W = fb2pp_system( w, D );
		Wf          = subproblem.td2fd( W );
	}

//	void test_gradient( const ivec& idxs, const rvec& coeffs, const rseq& h0, const rseq& f0 )
//	{
//		sbindexes = all_sbindexes.indirect( idxs );
//		sbvalues  = coeffs;
//		criterion::set_analysis_fb( h0 );
//		criterion::set_synthesis_fb( f0 );
//		cout << "Theoretical gradient = " << criterion::gradient() << endl;
//
//		int l  = h0.size();
//		int ht1 = h0.t1();
//		int ft1 = f0.t1();
//		auto fun = [&]( const rvec& x )
//		{
//			int n = idxs.size();
//			rvec coeffs = x.range( 0, n-1 );
//			cseq h0     = x.range( n, n+l-1 ); h0.t1(ht1);
//			cseq f0     = x.range( n+l, n+2*l-1 ); f0.t1(ft1);
//			logarithmic sbprob = *this;
//			sbprob.set_analysis_fb(h0);
//			sbprob.set_synthesis_fb(f0);
//			sbprob.sbvalues  = coeffs;
//			return sbprob.errfun();
//		};
//		cout << "Numerical gradient   = " << gradient(fun)(rvec{coeffs,h0,f0}) << endl;
//	}
};

/// @}
}
}
}
#endif
