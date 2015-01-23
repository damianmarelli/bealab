// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

#include <bealab/extensions/signal/sbsysapp.hpp>

namespace bealab
{
namespace signal
{
namespace sbsysapp
{
//==============================================================================
// Struct sbindex
//==============================================================================
bool sbindex::operator==( const sbindex& idx ) const
{
	if( idx.m == m && idx.n == n && idx.t == t && idx.real_f == real_f && idx.num_f == num_f )
		return true;
	else
		return false;
}

bool sbindex::operator!=( const sbindex& idx ) const
{
	return !this->operator==(idx);
}

bool sbindex::operator<( const sbindex& idx ) const
{
	if( m != idx.m )
		return m < idx.m;
	if( n != idx.n )
		return n < idx.n;
	if( t != idx.t )
		return t < idx.t;
	if( real_f != idx.real_f )
		return real_f < idx.real_f;
	if( num_f != idx.num_f )
		return num_f < idx.num_f;
	return false;
}

sbindex sbindex::conjugate( int M ) const
{
	int mc = mod( -m, M );
	int nc = mod( -n, M );
	return sbindex{ mc, nc, t, real_f, num_f };
}

bool sbindex::is_selfconjugate( int M ) const
{
	if( *this == this->conjugate(M) )
		return true;
	else
		return false;
}

//==============================================================================
// Class sbindex_pool
//==============================================================================
int sbindex_pool::sbi2int( const sbindex& idx )
{
	auto it = find( sbindexes.begin(), sbindexes.end(), idx );
	return it - sbindexes.begin();
}

void sbindex_pool::set_sbindexes( const Vec<sbindex>& idxs )
{
	sbindexes.resize(0);
	int I = idxs.size();
	for( int i = 0; i < I; i++ ) {

		// Check if the index is already listed
		sbindex idx = idxs(i);
		if( sbi2int(idx) < (int)sbindexes.size() ) {
//				cout << "Removed repeated index " << idx << endl;
			continue;
		}

		// In the case of real_target==true, check if the conjugate index is
		// already in
		sbindex idxc = idx.conjugate( M );
		if( real_target_f && sbi2int(idxc) < (int)sbindexes.size() ) {
//				cout << "Removed conjugate index " << idx << endl;
			continue;
		}

		// In the case of real_target==true, check if it is imaginary and from
		// subband 0 or M/2;
		if( real_target_f && !idx.real_f && idx.m == idx.n &&
				(idx.m == 0 || idx.m == M/2.) ) {
//				cout << "Removed imaginary chotix index " << idx << endl;
			continue;
		}

		sbindexes = Vec<sbindex>{ sbindexes, {idx} };
	}
}

void sbindex_pool::set_sbindexes( const Mat<cseq>& MA, const Mat<cseq>& AR )
{
	// Init
	assert( (int)MA.size1() == M && (int)MA.size2() == M );
	assert( (int)AR.size1() == M && (int)AR.size2() == M );
	Vec<sbindex> idxs;

	// Set numerator coefficients
	for( int m = 0; m < M; m++ )
		for( int n = 0; n < M; n++ ) {
			int T1 = MA(m,n).t1();
			int T2 = MA(m,n).t2();
			for( int t = T1; t <= T2; t++ )
				idxs = Vec<sbindex>{ idxs, { sbindex{m,n,t,true,true},
						sbindex{m,n,t,false,true} } };
		}

	// Set denominator coefficients
	for( int n = 0; n < M; n++ )
		for( int m = 0; m < M; m++ ) {
			int T1 = AR(m,n).t1();
			int T2 = AR(m,n).t2();
			for( int t = T1; t <= T2; t++ )
				idxs = Vec<sbindex>{ idxs, { sbindex{m,n,t,true,false},
						sbindex{m,n,t,false,false} } };
		}

	set_sbindexes( idxs );
}

void sbindex_pool::set_sbindexes( int noffdiags, int t1, int t2, bool den_f )
{
	assert( noffdiags >= 0 );
	Mat<cseq> MA(M,M), AR(M,M);
	int d1 = -min( noffdiags, M/2 - 1 );
	int d2 =  min( noffdiags, M/2 );
	for( int m = 0; m < M; m++ )
		for( int no = m-d1; no <= m+d2; no++ ){
			int n     = mod( no, M );
			MA(m,n) = rseq( ones(t2-t1+1), t1 );
			if(den_f)
				AR(m,n) = rseq( ones(t2), 1 );
		}
	set_sbindexes( MA, AR );
}

int sbindex_pool::sparsity() const
{
	int sp = 0;
	int N  = sbindexes.size();
	for( int n = 0; n < N; n++ ) {
		sp++;
		sbindex sbi = sbindexes(n);
		if( sbi != sbi.conjugate(M) )
			sp++;
	}
	return sp;
}

//==============================================================================
// Class sbmodel
//==============================================================================
int sbmodel::inv_ImP_irlength( const Vec<cmat>& P ) const
{
	int L = P.size() - 1;
	int N = L*M;
	cmat T = zeros(N,N);
	for( int l = 0; l < L; l++ )
		T( range(0,M), range(l*M,(l+1)*M) ) = P(l+1);
	for( int n = 0; n < N-M; n++ )
		T(n+M,n) = 1;
	auto DV = eig( T );
	cmat D  = get<0>( DV );
	double dommode = max(abs(diag(D)));
	if( dommode > 1 )
		error("irlength() - Unstable mode");
	return ceil( log(.001)/log(dommode) );
}

Mat<cseq> sbmodel::inv_ImP( const Mat<cseq>& P_ ) const
{
	// Form ARMA model
	Seq<cmat> I = Vec<cmat>{eye(M)};
	Seq<cmat> B = I;
	Seq<cmat> P = ms2sm(P_);
	if( P.size() == 0 )
		return sm2ms(I);
	Seq<cmat> A = I - P;
	control::arma<cmat> H( B.vec(), A.vec() );

	// Compute IR length
	int L = inv_ImP_irlength( P.vec() );

	// Impulse input
	Seq<cmat> X(L);
	X(0) = eye(M);
	for( int l = 1; l < L; l++ )
		X(l) = zeros(M,M);

	// Impulse response
	Seq<cmat> Ir = H(X);
	return sm2ms(Ir);

//		// Verify
//		Mat<cseq> Y_ = sm2ms(Ir);
//		Mat<cseq> I_(M,M);
//		for( int m = 0; m < M; m++ )
//			I_(m,m) = rseq{1};
//		cout << "inv_ImP() - error = " << norm( I_ - (I_-P_)*Y_ ) << endl;
//		return Y_;
}

sbmodel::models_t sbmodel::models() const
{
	// Build subband model
	models_t sbm;
	sbm.S.resize(M,M);
	sbm.P.resize(M,M);
	int I = sbindexes.size();
	for( int i = 0; i < I; i++ ) {
		sbindex si    = sbindexes(i);
		sbindex sic   = si.conjugate( M );
		complex c = si.real_f ? complex(1) : j;
		if(si.num_f) {
			sbm.S(si.m,si.n)(si.t) += c * sbvalues(i);
			if( real_target_f && si != sic )
				sbm.S(sic.m,sic.n)(si.t) += conj(c) * sbvalues(i);
		}
		else {
			sbm.P(si.m,si.n)(si.t) += c * sbvalues(i);
			if( real_target_f && si != sic )
				sbm.P(sic.m,sic.n)(si.t) += conj(c) * sbvalues(i);
		}
	}

	return sbm;
}

Mat<cseq> sbmodel::impulse_response() const
{
	models_t sbm = models();
	Mat<cseq> X  = inv_ImP( sbm.P );
	return X * sbm.S;
}

//==============================================================================
// Class fdcriterion
//==============================================================================
Vec<cvec> fdcriterion::get_atom_sbm_( const sbindex& sbi )
{
	// Form the polyphase representation of the atom
	Mat<cseq> E = outer_prod( FA.column(sbi.m), H.row(sbi.n) );
	E           = entrywise( [&sbi]( const cseq& x ){ return time_shift(x,sbi.t); } )( E );

	// Modify E if the target is real and the index is not self-conjugate
	sbindex sbic = sbi.conjugate( M );
	if( real_target_f && sbi != sbic ) {
		if(sbi.real_f)
			E =  2 * real( E );
		else
			E = -2 * imag( E );
	}

	// Return the atom
	return entrywise( dtft )( pp2fb_system( E ) );
}

void fdcriterion::reset_sbatoms()
{
	real_atoms_sbm.clear();
	imag_atoms_sbm.clear();
}

void fdcriterion::precompute_fbatoms()
{
	Mat<cseq> S = sbm.models().S;
	FASWA       = sm2ms( ms2sm( FA * S ) * Seq<cmat>(Vec<cmat>{M*idft(M)}) );
	WSH         = sm2ms( Seq<cmat>(Vec<cmat>{dft(M)}) * ms2sm( S * H ) );
}

Vec<cvec>& fdcriterion::get_atom_sbm( const sbindex& sbi )
{
	// Retrieve the fullband frequency response of the corresponding atom
	Vec<cvec>& e = sbi.real_f ? real_atoms_sbm(sbi.m,sbi.n)(sbi.t):
								imag_atoms_sbm(sbi.m,sbi.n)(sbi.t);

	// Fill the corresponding entry if it is empty
	if( sum(sum(abs(e))) == 0 )
		e = get_atom_sbm_( sbi );

	// Return the atom
	return e;
}

Vec<cvec> fdcriterion::get_atom_analysis( int idx )
{
	Mat<cseq> E(D,D);
	int i       = mod( idx, M );
	int j       = mod( idx, D );
	int d       = floor( (double)idx / D );
	E.column(j) = entrywise( [d](const cseq& x){ return time_shift(x,d); } )( Vec<cseq>(FASWA.column(i)) );
	return entrywise( dtft )( pp2fb_system( E ) );
}

Vec<cvec> fdcriterion::get_atom_synthesis( int idx )
{
	Mat<cseq> E(D,D);
	int i    = mod( idx, D );
	int j    = mod( idx, M );
	int d    = -floor( (double)idx / D );
	E.row(i) = entrywise( [d](const cseq& x){ return time_shift(x,d); } )( Vec<cseq>(WSH.row(j)) );
	return entrywise( dtft )( pp2fb_system( E ) );
}

fdcriterion::fdcriterion( int M, int D, int N ) :
	sbapprox(M,D),
	real_atoms_sbm(M,M), imag_atoms_sbm(M,M), Ndtft(N), wf2(ones(N))
{
	dtft = [this]( const cseq& x )
	{
		return bealab::dtft( x, this->Ndtft );
	};
}

void fdcriterion::set_target( const Mat<cseq>& G )
{
	sbapprox::set_target(G);
	gf = entrywise(dtft)( pp2fb_system(G) );
}

void fdcriterion::set_analysis_fb( const cseq& h0 )
{
	sbapprox::set_analysis_fb(real(h0));
	reset_sbatoms();
}

void fdcriterion::set_synthesis_fb( const cseq& f0 )
{
	sbapprox::set_synthesis_fb(real(f0));
	reset_sbatoms();
}

void fdcriterion::set_spectral_weight( const rseq& w_ )
{
	w   = w_;
	wf2 = real( pow( abs( dtft(w) ), 2 ) );
}

//==============================================================================
// Class logarithmic
//==============================================================================
rvec logarithmic::gradient_aux( const Vec<Vec<cvec>>& atoms )
{
	// Frequency response of the approximation
	Vec<cseq> gh  = impulse_responses();
	Vec<cvec> ghf = entrywise(dtft)( gh );

	// Pre-compute X = conj(\tilde{c}) / \hat{g}
	Vec<cvec> X(D);
	for( int d = 0; d < D; d++ ) {
		cvec ctfd = logspec( element_div( gf(d), ghf(d) ) );
		X(d) = element_div( conj(ctfd), ghf(d) );
	}

	// Compute the gradient entry-by-entry
	double K = -40/(Ndtft*D*log(10));
	int I    = atoms.size();
	rvec grad(I);
	for( int i = 0; i < I; i++ ) {

		// Get the i-th atom
		const Vec<cvec>& ei = atoms(i);

		// Gradient entry
		double gi = 0;
		for( int d = 0; d < D; d++ ) {
			rvec t1 = element_prod( real(X(d)), real(ei(d)) );
			rvec t2 = element_prod( imag(X(d)), imag(ei(d)) );
			rvec t3 = t1 - t2;
			rvec t4 = element_prod( wf2, t3 );
			gi     += sum( t4 );
		}

		// Set the gradient entry
		grad(i) = K * gi;
	}

	// Return gradient
	return grad;
}

logarithmic::logarithmic( int M, int D, int N ) :
	sbapprox(M,D),
	fdcriterion(M,D,N),
	fig_amp("Spectra"), fig_phase("Phase")
{
	logspec = []( const cvec& x ) -> cvec
	{
		return 20*log10( x );
	};
}

double logarithmic::errfun()
{
	Mat<cseq> Gh  = polyphase();
	Vec<cseq> gh  = pp2fb_system( Gh );
	Vec<cvec> ghf = entrywise(dtft)( gh );
	Vec<cvec> ctf(D);
	for( int d = 0; d < D; d++ )
		ctf(d) = logspec( element_div( gf(d), ghf(d) ) );
	rvec x1 = sum( real( pow( abs(ctf), 2 ) ) );
	rvec x2 = element_prod( wf2, x1 );
	return 1./(Ndtft*D) * sum( x2 );
}

rvec logarithmic::gradient_sbm()
{
	// Form the vector of atoms
	int I = sbm.sbindexes.size();
	Vec<Vec<cvec>> atoms(I);
	for( int i = 0; i < I; i++ )
		atoms(i) = get_atom_sbm( sbm.sbindexes(i) );

	// Gradient
	rvec grad = gradient_aux( atoms );

	// Return the gradient
	return grad;
}

logarithmic::gradient_t logarithmic::gradient()
{
	// Current approximation
	precompute_fbatoms();

	// Gradient
	gradient_t grad;

	// Gradient with respect to the subband model
	grad.sbm = gradient_sbm();

	// Gradient with respect to the analysis window
	{
		// Form the vector of analysis atoms
		int I = h0().size();
		Vec<Vec<cvec>> atoms(I);
		for( int i = 0; i < I; i++ )
			atoms(i) = get_atom_analysis( h0().t1()+i );
		grad.h0 = gradient_aux( atoms );
	}

	// Gradient with respect to the synthesis window
	{
		// Form the vector of analysis atoms
		int I = f0().size();
		Vec<Vec<cvec>> atoms(I);
		for( int i = 0; i < I; i++ )
			atoms(i) = get_atom_synthesis( f0().t1()+i );
		grad.f0 = gradient_aux( atoms );
	}

	// Return the gradient
	return grad;
}

void logarithmic::plotfun()
{
	// Functors
	int N  = 2 * Ndtft;
	rvec w = linspace( -pi, pi, N+1 );
	w      = w( range( 0, w.size()-1 ) );
	int i0 = find( w.begin(), w.end(), 0 ) - w.begin();
	function<rvec(const rvec&)> centerphase = [i0](const rvec& x){ return x - x(i0)*ones(x.size()); };
	function<cvec(const cvec&)> ftshift     = [](const cvec& x){ return bealab::ftshift(x); };
	function<cvec(const cvec&)> dbspec      = [](const cvec& x ){ return 20*log10(x); };
	function<cvec(const cseq&)> dtft        = [N](const cseq& x ){ return bealab::dtft(x,N); };

	// Target frequency response
	Vec<cseq> g   = pp2fb_system( G );
	Vec<cvec> gf  = entrywise(dtft)( g );
	Vec<cvec> cf  = entrywise(ftshift)( entrywise(dbspec)( gf ) );
	Vec<rvec> cfr = real(cf);
	Vec<rvec> cfi = real(-i*cf);
	cfi = 180/pi * entrywise(unwrap)( log(10)/20 * cfi );
	cfi = entrywise(centerphase)( cfi );

	// Approximation frequency response
	Mat<cseq> Gh   = polyphase();
	Vec<cseq> gh   = pp2fb_system( Gh );
	Vec<cvec> ghf  = entrywise(dtft)( gh );
	Vec<cvec> chf  = entrywise(ftshift)( entrywise(dbspec)(ghf) );
	Vec<rvec> chfr = real(chf);
	Vec<rvec> chfi = real(-i*chf);
	chfi = 180/pi * entrywise(unwrap)( log(10)/20 * chfi );
	chfi = entrywise(centerphase)( chfi );

//		// Error
//		Vec<cvec> gtf(D);
//		for( int d = 0; d < D; d++ )
//			gtf(d) = element_div( gf(d), ghf(d) );
//		Vec<cvec> ctf = gtf.apply(logspec).apply( ftshift );

	// Plot
	fig_amp.clear().overlap().rangex(-pi,pi).
			plot( w, abs(cfr-chfr),   grey ).
			plot( w, chfr,  blue ).
			plot( w, cfr,   red );
	fig_phase.clear().overlap().rangex(-pi,pi).
			plot( w, abs(cfi-chfi),   grey ).
			plot( w, chfi,  blue ).
			plot( w, cfi,   red );
}

//==============================================================================
// Class logarithmic_nophase
//==============================================================================
rvec logarithmic_nophase::gradient_aux( const Vec<Vec<cvec>>& atoms )
{
	// Frequency response of the approximation
	Vec<cseq> gh  = impulse_responses();
	Vec<cvec> ghf = entrywise(dtft)( gh );

	// Some pre-computing
	Vec<rvec> C(D);
	for( int d = 0; d < D; d++ ) {
		rvec t1   = real( logspec( element_div( gf(d), ghf(d) ) ) );
		rvec t2   = element_prod( wf2, t1 );
		rvec t3   = real( pow( abs( ghf(d) ), 2 ) );
		C(d)      = element_div( t2, t3 );
	}

	// Compute the gradient entry-by-entry
	double K = -40/(Ndtft*D*log(10));
	int I    = atoms.size();
	rvec grad(I);
	for( int i = 0; i < I; i++ ) {

		// Get the i-th atom
		const Vec<cvec>& ei = atoms(i);

		// Gradient entry
		double gi = 0;
		for( int d = 0; d < D; d++ ) {
			rvec t1 = element_prod( real(ghf(d)), real(ei(d)) );
			rvec t2 = element_prod( imag(ghf(d)), imag(ei(d)) );
			rvec t3 = t1 + t2;
			rvec t4 = element_prod( C(d), t3 );
			gi     += sum( t4 );
		}

		// Set the gradient entry
		grad(i) = K * gi;
	}

	// Return gradient
	return grad;
}

logarithmic_nophase::logarithmic_nophase( int M, int D, int N ) :
	sbapprox(M,D),
	logarithmic(M,D,N)
{
	logspec = []( const cvec& x ) -> cvec
	{
		return 20*log10( abs( x ) );
	};
}

//==============================================================================
// Class greedy
//==============================================================================
int greedy::number_of_coefficients( const Mat<cseq>& S )
{
	return norm(real(S),0) + norm(imag(S),0);
}

bool greedy::forcestop()
{
	if( (int)sbm.sparsity() >= sparse::greedy<double>::max_sparsity-1 )
		return true;
	if( sbm.sbindexes.size() >= sbindexes.size() )
		return true;
	return false;
}

void greedy::plotfun()
{
	// Plot from the base
	sparse::greedy<double>::plotfun();

	// Functor to parse a coefficient
	rvec mr, nr, tr, mi, ni, ti;
	auto parse = [&]( const sbindex& idx )
	{
		if( idx.real_f ) {
			mr = { mr, {(double)idx.m} };
			nr = { nr, {(double)idx.n} };
			tr = { tr, {(double)idx.t} };
		}
		else{
			mi = { mi, {(double)idx.m} };
			ni = { ni, {(double)idx.n} };
			ti = { ti, {(double)idx.t} };
		}
	};

	// Plot the coefficients
	int I = this->indexes.size();
	for( int i = 0; i < I; i++ ) {
		sbindex idx  = sbindexes( this->indexes(i) );
		sbindex idxc = idx.conjugate( M );
		parse( idx );
		if( real_target_f && idx != idxc )
			parse( idxc );
	}
	fig.clear().overlap()
			.plot( mr, nr, tr, star, noline, blue, "Real" )
			.plot( mi, ni, ti, circle, noline, red, "Imaginary" );
}

greedy::greedy( int M, int D ) :
	sbapprox(M,D),
	fig("coefficients") {}

void greedy::approximate()
{
	sparse::greedy<double>::approximate();
	sbm.sbindexes = sbindexes(indirect(this->indexes));
	sbm.sbvalues  = this->coefficients;
}

//==============================================================================
// Namespace kronecker
//==============================================================================
namespace kronecker
{
//==============================================================================
// Class linear
//==============================================================================
void linear::plotfun( const ivec& idxs, const rvec& coeffs )
{
	// Plotting tools
	int N  = 1024;
	rvec w = linspace( -pi, pi, N+1 );
	w      = w( range( 0, w.size()-1 ) );
	function<rvec(const rvec&)> ftshift = [](const rvec& x ){ return bealab::ftshift(x); };
	function<rvec(const cseq&)> spec    = [N](const cseq& x ) { return abs( dtft(x,N) ); };

	// Target frequency response
	rvec gw = ftshift( spec( pp2fb_system(G)(0) ) );

	// Approximation frequency response
	sbm.sbindexes = sbindexes(indirect(idxs));
	sbm.sbvalues  = coeffs;
	Vec<rseq> gh  = real(impulse_responses());
	Vec<rvec> ghw = entrywise(ftshift)( entrywise(spec)(gh) );

	// Plot
	fig.clear().overlap().rangex(-pi,pi).
			plot( w, gw ).
			plot( w, ghw, blue );
}

//==============================================================================
// Class time_domain
//==============================================================================
Mat<cseq> time_domain::synthesis_( const ivec& idxs, const rvec& coeffs )
{
	sbm.sbindexes = sbindexes(indirect(idxs));
	sbm.sbvalues  = coeffs;
	return polyphase();
}

rvec time_domain::analysis_( const Mat<cseq>& G )
{
//		if( analysis_1by1_flag )
//			return analysis_1by1( G );
//		else
//			return analysis_whole( G );
	const Mat<cseq>& S = F * G * adjoint(H);
	int I = sbindexes.size();
	rvec x(I);
	for( int i = 0; i < I; i++ ) {
		sbindex si = sbindexes(i);
		if(si.real_f)
			x(i) = real( S(si.m,si.n)(si.t) );
		else
			x(i) = imag( S(si.m,si.n)(si.t) );

		// Correction for real target
		if( real_target_f && !si.is_selfconjugate(M) )
			x(i) *= 2;
	}
	return x;
}

double time_domain::inner_product_atoms_( int i, int j )
{
	sbindex sbii = sbindexes(i);
	sbindex sbij = sbindexes(j);
	complex s = (FF(sbij.m,sbii.m) * HH(sbii.n,sbij.n))(sbij.t-sbii.t);
	double ip;
	if( sbii.real_f && sbij.real_f )
		ip = real(s);
	else if( sbii.real_f && !sbij.real_f )
		ip = imag(s);
	else if( !sbii.real_f && sbij.real_f )
		ip = -imag(s);
	else
		ip = real(s);

	// Correction for real target
	if( real_target_f &&
			( !sbii.is_selfconjugate(M) || !sbij.is_selfconjugate(M) ) )
		ip *= 2;

	return ip;
}

ivec time_domain::next_indexes()
{
	ivec idxs = orthogonal_matching_pursuit::next_indexes();
	if(this->trace)
		cout << "Chosen index: " << sbindexes( idxs(idxs.size()-1) ) << endl;
	return idxs;
}

void time_domain::plotfun()
{
	greedy::plotfun();
	linear::plotfun( indexes, coefficients );
}

time_domain::time_domain( int M, int D ) :
	sbapprox(M,D),
	greedy(M,D),
	linear(M,D)
{
	// Define the functors of orthogonal_matching_pursuit
	synthesis = [this]( const ivec& idxs, const rvec& coeffs )
	{
		return this->synthesis_( idxs, coeffs );
	};

	analysis = [this]( const Mat<cseq>& G )
	{
		return this->analysis_( G );
	};

	inner_product_atoms = [this]( int i, int j )
	{
		return this->inner_product_atoms_( i, j );
	};
}

void time_domain::set_target( const Mat<cseq>& G )
{
	sbapprox::set_target( G );
	target = G;
}

void time_domain::set_analysis_fb( const Mat<cseq>& H )
{
	sbapprox::set_analysis_fb( H );
	HH = H * adjoint(H);
}

void time_domain::set_synthesis_fb( const Mat<cseq>& F )
{
	sbapprox::set_synthesis_fb( F );
	FF = F * adjoint(F);
}

void time_domain::set_sbindexes( const Vec<sbindex>& idxs )
{
	Vec<sbindex> idxs_num;
	int I = idxs.size();
	for( int i = 0; i < I; i++ )
		if(idxs(i).num_f)
			idxs_num = { idxs_num, {idxs(i)} };
	greedy::set_sbindexes( idxs_num );
}

void time_domain::approximate()
{
	orthogonal_matching_pursuit::initialization();
	greedy::approximate();
}

//==============================================================================
// Class frequency_domain
//==============================================================================
Vec<cmat> frequency_domain::synthesis_( const ivec& idxs, const rvec& coeffs )
{
	sbm.sbindexes = sbindexes(indirect(idxs));
	sbm.sbvalues  = coeffs;
	Mat<cseq> S   = sbm.models().S;
	Vec<cmat> Sf  = td2fd( S );
	Vec<cmat> Ghf = element_prod( element_prod( adjoint(Ff), Sf ), Hf );
	return Ghf;
}

rvec frequency_domain::analysis_( const Vec<cmat>& Gf )
{
	// Compute the Mat<cseq> of complex coefficients
	Vec<cmat> Sf  = element_prod( element_prod( Ff, Gf ), adjoint(Hf) );
	Mat<cseq> S   = fd2td( Sf );

	// Compute the inner products
	int I = sbindexes.size();
	rvec x(I);
	for( int i = 0; i < I; i++ ) {
		sbindex si = sbindexes(i);
		if(si.real_f)
			x(i) = real( S(si.m,si.n)(si.t) );
		else
			x(i) = imag( S(si.m,si.n)(si.t) );

		// Correction for real target
		if( real_target_f && !si.is_selfconjugate(M) )
			x(i) *= 2;
	}
	return x;
}

double frequency_domain::inner_product_atoms_( int i, int j )
{
	sbindex sbii  = sbindexes(i);
	sbindex sbij  = sbindexes(j);
	cvec u    = element_prod( FFf(sbij.m,sbii.m), HHf(sbii.n,sbij.n) );
	cseq v    = idtft( u );
	complex s = v(sbij.t-sbii.t);
	double ip;
	if( sbii.real_f && sbij.real_f )
		ip = real(s);
	else if( sbii.real_f && !sbij.real_f )
		ip = imag(s);
	else if( !sbii.real_f && sbij.real_f )
		ip = -imag(s);
	else
		ip = real(s);

	// Correction for real target
	if( real_target_f &&
			( !sbii.is_selfconjugate(M) || !sbij.is_selfconjugate(M) ) )
		ip *= 2;

	return ip;
}

ivec frequency_domain::next_indexes()
{
	ivec idxs = orthogonal_matching_pursuit::next_indexes();
	if(this->trace)
		cout << "Chosen index: " << sbindexes(idxs(idxs.size()-1)) << endl;
	return idxs;
}

/// Virtual function override
void frequency_domain::plotfun()
{
	greedy::plotfun();
	linear::plotfun( indexes, coefficients );
}

frequency_domain::frequency_domain( int M, int D, int N ) :
	sbapprox(M,D),
	greedy(M,D),
	linear(M,D)
{
	// Initialize functors
	dtft  = [N]( const cseq& x ) { return bealab::dtft(x,N); };
	idtft = []( const cvec& x ) { return bealab::idtft(x); };

	td2fd = [this]( const Mat<cseq>& G )
	{
		auto tmp = entrywise(this->dtft)( G );
		return ms2sm( Mat<cseq>( tmp ) ).vec();
	};

	fd2td = [this]( const Vec<cmat>& Gf )
	{
		auto tmp = entrywise(this->idtft)( vm2mv(Gf) );
		return tmp;
	};

	// Functors of orthogonal_matching_pursuit
	synthesis = [this]( const ivec& idxs, const rvec& coeffs )
	{
		return this->synthesis_( idxs, coeffs );
	};

	analysis = [this]( const Vec<cmat>& Gf )
	{
		return this->analysis_( Gf );
	};

	inner_product_atoms = [this]( int i, int j )
	{
		return this->inner_product_atoms_( i, j );
	};
}

void frequency_domain::set_target( const Vec<cmat>& Gf )
{
	target = Gf;
}

void frequency_domain::set_target( const Mat<cseq>& G )
{
	sbapprox::set_target( G );
	Vec<cmat> Gf = td2fd( G );
	set_target( Gf );
}

void frequency_domain::set_analysis_fb( const Vec<cmat>& Hf_ )
{
	Hf             = Hf_;
	Vec<cmat> HHf1 = element_prod( Hf, adjoint(Hf) );
	HHf            = vm2mv( HHf1 );
}

void frequency_domain::set_analysis_fb( const Mat<cseq>& H )
{
	sbapprox::set_analysis_fb( H );
	Vec<cmat> Hf_  = td2fd( H );
	set_analysis_fb( Hf_ );
}

void frequency_domain::set_synthesis_fb( const Vec<cmat>& Ff_ )
{
	Ff             = Ff_;
	Vec<cmat> FFf1 = element_prod( Ff, adjoint(Ff) );
	FFf            = vm2mv( FFf1 );
}

void frequency_domain::set_synthesis_fb( const Mat<cseq>& F )
{
	sbapprox::set_synthesis_fb( F );
	Vec<cmat> Ff_  = td2fd( F );
	set_synthesis_fb( Ff_ );
}

void frequency_domain::set_sbindexes( const Vec<sbindex>& idxs )
{
	Vec<sbindex> idxs_num;
	int I = idxs.size();
	for( int i = 0; i < I; i++ )
		if(idxs(i).num_f)
			idxs_num = { idxs_num, {idxs(i)} };
	greedy::set_sbindexes( idxs_num );
}

void frequency_domain::approximate()
{
	orthogonal_matching_pursuit::initialization();
	greedy::approximate();
}

}
}
}
}
