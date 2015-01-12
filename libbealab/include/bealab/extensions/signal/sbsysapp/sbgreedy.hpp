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

	/// Count the number of non-zero coefficients in a Mat<cseq>
	int number_of_coefficients( const Mat<cseq>& S );

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
	public sparse::orthogonal_matching_pursuit<double,Mat<cseq>> {

private:

	Mat<cseq> HH;																///< Precomputed H*H.star()
	Mat<cseq> FF;																///< Precomputed F*F.star()

	/// Callback function
	Mat<cseq> synthesis_( const ivec& idxs, const rvec& coeffs );

	/// Callback function
	rvec analysis_( const Mat<cseq>& G );

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
	void set_target( const Mat<cseq>& G ) override;

	/// Virtual function override
	void set_analysis_fb( const Mat<cseq>& H ) override;

	/// Virtual function override
	void set_synthesis_fb( const Mat<cseq>& F ) override;

	/// Virtual function override
	void set_sbindexes( const Vec<sbindex>& idxs ) override;

	/// Virtual function override
	void approximate() override;
};

/// Frequency domain.
class frequency_domain :
	public greedy,
	public linear,
	public sparse::orthogonal_matching_pursuit<double,Vec<cmat>> {

	Mat<cvec> HHf;																///< Precomputed H*H.star()
	Mat<cvec> FFf;																///< Precomputed F*F.star()
	Vec<cmat> Hf;																///< Frequency version of H
	Vec<cmat> Ff;																///< Frequency version of F

	/// Callback function
	Vec<cmat> synthesis_( const ivec& idxs, const rvec& coeffs );

	/// Callback function
	rvec analysis_( const Vec<cmat>& Gf );

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

	function<Vec<cmat>(const Mat<cseq>&)> td2fd;
	function<Mat<cseq>(const Vec<cmat>&)> fd2td;

	// Imports from subband
	using sbapprox::set_target;
	using sbapprox::set_analysis_fb;
	using sbapprox::set_synthesis_fb;
	using sbapprox::set_sbindexes;

	/// Constructor
	frequency_domain( int M, int D, int N );

	/// Set the target
	void set_target( const Vec<cmat>& Gf );

	/// Virtual function override
	void set_target( const Mat<cseq>& G ) override;

	/// Set the analysis filterbank
	void set_analysis_fb( const Vec<cmat>& Hf_ );

	/// Virtual function override
	void set_analysis_fb( const Mat<cseq>& H ) override;

	/// Set the synthesis filterbank
	void set_synthesis_fb( const Vec<cmat>& Ff_ );

	/// Virtual function override
	void set_synthesis_fb( const Mat<cseq>& F ) override;

	/// Virtual function override
	void set_sbindexes( const Vec<sbindex>& idxs ) override;

	/// Virtual function override
	void approximate() override;
};

}

/// Subband system approximation using orthogonal matching pursuit
template<class criterion>
class ompursuit : public kronecker::time_domain {

	criterion crit;
	rseq w;																		///< Spectral weight
	Mat<rseq> W;																///< Polyphase representation of the spectral weight factor

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
//	void set_target( const Mat<cseq>& G ) override
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
		Mat<cseq> G_ = G;
		Mat<cseq> H_ = H;
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
	Mat<rseq> W;																///< Polyphase representation of the spectral weight factor
	Vec<cmat> Rf;																///< Frequency-polyphase representation of the modified spectral weights

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
		Mat<cseq> Gh  = polyphase();
		Mat<cseq> Gt  = G - Gh;
		Vec<cmat> Gtf = td2fd( Gt );

		// Compute log difference
		Vec<cvec> gf  = entrywise(dtft)( pp2fb_system(G) );
		Vec<cvec> ghf = entrywise(dtft)( pp2fb_system( Gh ) );
		Vec<cvec> ctf(D);
		for( int d = 0; d < D; d++ )
			ctf(d) = crit.logspec( element_div( crit.gf(d), ghf(d) ) );
		Vec<cseq> ct   = entrywise( [this](const cvec& x){ return idtft(x); } )( ctf );
		Mat<cseq> Ct   = fb2pp_system( ct );
		Vec<cmat> Ctf  = td2fd( Ct );
		Vec<cmat> Wf   = td2fd( W );
		Vec<cmat> CtWf = element_prod( Ctf, Wf );

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
		Mat<cseq> G_ = G;
		Mat<cseq> H_ = H;

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
			Vec<cmat> Hf  = td2fd( H );
			Vec<cmat> Gf  = td2fd( G );
			Vec<cmat> Gzf = element_prod( Gf, Rf );
			Vec<cmat> Hzf = element_prod( Hf, Rf );
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
	public sparse::gradient_pursuit<Mat<cseq>> {

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
	Vec<cmat> Gf;																///< Frequency-polyphase representation of the target
	Vec<cmat> Hf;																///< Frequency-polyphase representation of the analysis filterbank
	Vec<cmat> Wf;																///< Frequency-polyphase representation of the nominal spectral weight
	Vec<cmat> Rf;																///< Frequency-polyphase representation of the modified spectral weights
	bool gradientpursuit_f;														///< Flag to indicate if the next indexes are chosen using gradient pursuit

	/// Update the auxiliar variables Zf and CtWf
	void update_auxiliar_variables()
	{
//		// Compute linear difference
//		Vec<cseq> gh  = pp2fb_system( approximation );
//		Vec<cvec> ghf = gh.apply( dtft );
//		Vec<cvec> gtf = gf - ghf;
//
//		// Compute log difference
//		Vec<cvec> ctf(D);
//		for( int d = 0; d < D; d++ )
//			ctf(d) = logspec( element_div( gf(d), ghf(d) ) );
//
//		// Spectral compensation weights
//		Vec<cvec> af(D);
//		for( int d = 0; d < D; d++ )
//			af(d) = element_div( ctf(d), gtf(d) );
//
//		// Average weight
//		rvec af1 = sum( abs( af ) ) / D;
//		rvec zf1 = element_prod( af1, sqrt(wf) );
//		Vec<rvec> zf(D);
//		for( int d = 0; d < D; d++ )
//			zf(d) = zf1;
//		cseq z1      = idtft( zf1 );
//		Mat<cseq> Z  = fb2pp_system( z1, D );
//		Zf           = sp_indexes.td2fd( Z );
//
//		// Compute the spectral compensation for logarithmic amplitude
//		if( indexes.size() == 0 )
//			Zf  = Wf;

		// Compute linear difference
		Mat<cseq> Gh  = polyphase();
		Mat<cseq> Gt  = G - Gh;
		Vec<cmat> Gtf = subproblem.td2fd( Gt );

		// Compute log difference
		Vec<cseq> gh  = pp2fb_system( Gh );
		Vec<cvec> ghf = entrywise(this->dtft)( gh );
		Vec<cvec> ctf(D);
		for( int d = 0; d < D; d++ )
			ctf(d) = this->logspec( element_div( this->gf(d), ghf(d) ) );
		Vec<cseq> ct   = entrywise( [](const cvec& x){ return idtft(x); } )( ctf );
		Mat<cseq> Ct   = fb2pp_system( ct );
		Vec<cmat> Ctf  = subproblem.td2fd( Ct );
		Vec<cmat> CtWf = element_prod( Ctf, Wf );

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
		Vec<cmat> Gzf = element_prod( Gf, Rf );
		Vec<cmat> Hzf = element_prod( Hf, Rf );
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
		Mat<cseq> S   = sbm.models().S;
		for( int m = 0; m < M; m++ )
			if( S(m,m).size() == 0 )
				return false;
		return true;
	}

	/// Virtual function override
	void set_analysis_fb( const Mat<cseq>& H ) override
	{
		// Set the base
		greedy::set_analysis_fb(H);

		// Frequency-polyphase representation of the analysis filterbank
		Hf = subproblem.td2fd( H );
	}

	/// Virtual function override
	void set_synthesis_fb( const Mat<cseq>& F ) override
	{
		// Set the base
		greedy::set_synthesis_fb(F);

		// Compute frequency-polypase of the dual
		Vec<cmat> Ff = subproblem.td2fd( F );
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

//		Mat<cseq> R   = subproblem.fd2td( Rf );
//		Vec<cseq> r   = pp2fb_system( R );
//		Vec<cvec> rf  = entrywise(this->dtft)( r );
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
	void set_target( const Mat<cseq>& G ) override
	{
		// Set the bases
		criterion::set_target( G );

		// Set local variables
		Gf = subproblem.td2fd( G );
	}

	/// Virtual function override
	void set_sbindexes( const Vec<sbindex>& idxs ) override
	{
		// Set the bases
		greedy::set_sbindexes(idxs);
		gradient_pursuit::set_size(sbindexes.size());

		// Set the subproblem
		int N = idxs.size();
		Vec<sbindex> sidxs(N);
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
		Mat<cseq> W = fb2pp_system( w, D );
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
