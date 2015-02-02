/*******************************************************************************
 * This software is licensed under the BSD 3-Clause License with the possibility to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
 * You may not use this work except in compliance with the License.
 * You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
 * Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the file License for the specific language governing permissions and limitations under the License. 
 * If you wish to obtain a commercial license, please contact the authors via e-mail.
 *
 * Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)
 *******************************************************************************/
#include <bealab/scilib/stats.hpp>
#include <bealab/scilib/calculus.hpp>
#include <bealab/scilib/clustering.hpp>
#include <bealab/scilib/rootfind.hpp>
#include <bealab/scilib/arrays.hpp>
#include <bealab/extensions/communications/quantization.hpp>

namespace bealab
{
namespace comms
{

//============================================================================
// Quantizer class
//============================================================================
quantizer::quantizer( const rvec &dict_ ) : dict(dict_)
{
	bnd = qboundaries(dict);
}

int quantizer::encode( double y ) const
{
	return find(bnd < y).size()-1;
}

iseq quantizer::encode( const rseq &y ) const
{
	iseq i( y.size(), y.t1() );
	for( int t = i.t1(); t <= i.t2(); t++ )
		i(t) = encode( y(t) );
	return i;
}

double quantizer::decode( int i ) const
{
	return dict(i);
}

rseq quantizer::decode( const iseq &i ) const
{
	rseq z( i.size(), i.t1() );
	for( int t = i.t1(); t <= i.t2(); t++ )
		z(t) = decode( i(t) );
	return z;
}

double quantizer::operator()( double x ) const
{
	return decode( encode(x) );
}

rseq quantizer::operator()( const rseq &x ) const
{
	return decode( encode(x) );
}

//============================================================================
// Vector quantizer class
//============================================================================
vector_quantizer::vector_quantizer( const vec<rvec> &database, int dictsize )
{
	train( database, dictsize );
}

void vector_quantizer::train( const vec<rvec> &database, int dictsize )
{
	dict = kmeans( database, dictsize );
}

int vector_quantizer::encode( const rvec &x ) const
{
	int D = dict.size();
	rvec dist(D);
	for( int d = 0; d < D; d++ )
		dist(d) = norm( x - dict(d) );
	return min_index( dist );
}

iseq vector_quantizer::encode( const vec<rseq> &x ) const
{
	return encode( vs2sv(x) );
}

iseq vector_quantizer::encode( const sequence<rvec> &x ) const
{
	iseq y( x.size(), x.t1() );
	for( int t = x.t1(); t <= x.t2(); t++ )
		y(t) = encode(x(t));
	return y;
}

rvec vector_quantizer::decode( int i ) const
{
	return dict(i);
}

sequence<rvec> vector_quantizer::decode( const iseq &i ) const
{
	sequence<rvec> z( i.size(), i.t1() );
	for( int t = i.t1(); t <= i.t2(); t++ )
		z(t) = decode( i(t) );
	return z;
}

rvec vector_quantizer::operator()( const rvec &x ) const
{
	return decode( encode(x) );
}

vec<rseq> vector_quantizer::operator()( const vec<rseq> &x ) const
{
	return sv2vs( decode( encode(x) ) );
}

sequence<rvec> vector_quantizer::operator()( const sequence<rvec> &x ) const
{
	return decode( encode(x) );
}

//============================================================================
// Quantizer design
//============================================================================
rvec qboundaries( const rvec &dict )
{
	int		N   = dict.size();
	rvec	bnd(N+1);
	bnd(0) = -inf;
	for( int n = 1; n < N; n++ ) {
		double	m2 = dict(n);
		double	m1 = dict(n-1);
		bnd(n) = (m2 + m1) / 2;
	}
	bnd(N) = inf;

	return bnd;
}

rvec qdict( const function<double(double)>& dist, const rvec &bnd )
{
	int		N   = bnd.size() - 1;
	rvec	dict(N);
	for( int n = 0; n < N; n++ ) {
		double b2 = bnd(n+1);
		double b1 = bnd(n);
		function<double(double)> fun = [&dist](double x){return x*dist(x);};
		dict(n)   = integral( fun, b1, b2 ) / integral( dist, b1, b2 );
	}

	return dict;
}

/*
 * Computes a Lloyd's quantizer for an arbitrary distribution
 */
quantizer lloyd( const function<double(double)>& dist, const rvec &bnd0 )
{
	rvec bnd    = bnd0;
	rvec dict   = qdict( dist, bnd );
	double	err = qerror( dist, quantizer(dict) );
	double err0;
	do {
		bnd  = qboundaries( dict );
		dict = qdict( dist, bnd );
		err0 = err;
		err  = qerror( dist, quantizer(dict) );
	} while( (err0-err)/err0 > 1e-6 );

	quantizer qnt;
	qnt.bnd  = bnd;
	qnt.dict = dict;
	return qnt;
}

/*
 * Computes a Lloyd's quantizer for a normal distribution
 */
quantizer lloyd( double mean, double variance, int L )
{
	double M = 4*sqrt(variance);
	rvec bnd = linspace( mean-M, mean+M, L+1 );
	std::function<double(double)> dist = [mean, variance](double x)
	{
		return 1./sqrt(2*pi*variance) * exp( -pow(x-mean,2)/(2*variance) );
	};

	return lloyd( dist, bnd );
}

/*
 * Waterfilling to allocate bits for quantizing a vector
 */
rvec waterfilling( const rvec &s2, int B )
{
	int M = s2.size();
	rvec R = zeros(M);
	for( int i = 0; i < B; i++ ) {
		rvec D = element_prod( 16./3 * real(pow( 2, -2*R )), s2 );
		int idx = max_index( D );
		R(idx) += 1;
	}

	return R;
}

/*
 * Filterbank for KL decomposition
 */
vec<rseq> KL_filterbank( const rseq &rx, int N, rvec *power, double tol )
{
	// input eigenvectors
	rmat R = zeros(N,N);
	for( int i = 0; i < N; i++ )
		for( int j = 0; j < N; j++ )
			R(i,j) = rx(j-i);
	auto DV = eig( R );
	cmat D  = get<0>( DV );
	cmat V  = get<1>( DV );
	cvec d  = diag( D );

	// selection of number of eigenvalues
	int Nc = sum( static_cast<ivec>( abs(d) > tol * abs(d(0)) ) );

	// Analysis filterbank
	vec<rseq> h(Nc);
	for( int i = 0; i < Nc; i++ )
	    h(i) = rseq( real(flip(V.column(i))), 0 );

    if( power!=NULL ) {
    	power->resize( Nc );
    	for( int i = 0; i < Nc; i++ )
    		(*power)(i) =  abs(d(i));
    }

    return h;
}

//============================================================================
// Analysis of general distributions
//============================================================================
/*
 * Computes the quantization error for an arbitrary distribution and quantizer
 */
double qerror( const function<double(double)>& dist, const quantizer &qnt )
{
	double err = 0;
	int I = qnt.dict.size();
	for( int i = 0; i < I; i++ ) {
		double xo = qnt.dict(i);
		function<double(double)> fun = [&dist,xo](double x){return std::pow(x-xo,2)*dist(x);};
		err += integral( fun, qnt.bnd(i), qnt.bnd(i+1) );
	}
	return err;
}

/*
 * Computes \int_{b1}^{b2} dist(x) dx
 */
double I0( const function<double(double)>& dist, double b1, double b2 )
{
	return integral( dist, b1, b2 );
}

/*
 * Computes \int_{b1}^{b2} x*dist(x) dx
 */
double I1( const function<double(double)>& dist, double b1, double b2 )
{
	function<double(double)> fun = [&dist](double x){return x*dist(x);};
	return integral( fun, b1, b2 );
}

/*
 * Computes the conditional mean inside a quantization interval
 */
double cmean( const function<double(double)>& dist, double b1, double b2 )
{
	double	i1 = I1( dist, b1, b2 );
	double	i0 = I0( dist, b1, b2 );
	double	rv;
	if( i0 != 0 )
		rv = i1/i0;
	else if( b1 == -inf && b2 == inf )
		rv = 0;
	else if( b1 == -inf )
		rv = b2;
	else if( b2 == inf )
		rv = b1;
	else
		rv = (b1+b2)/2;

	return rv;
}

//============================================================================
// Analysis of Gaussian variables
//============================================================================
/*
 * Computes the power of a quantized signal
 */
double qpower( double mean, double var, const quantizer &qnt )
{
	double P = 0;
	int I = qnt.bnd.size() - 1;
	for( int i = 0; i < I; i++ ) {
		double	Pxo = I0( mean, var, qnt.bnd(i), qnt.bnd(i+1) );
		double	xo  = qnt.dict(i);
		P += std::pow(xo,2) * Pxo;
	}
	return P;
}

/*
 * Computes \int_{b1}^{b2} dist(x) dx
 */
double I1( double mean, double var, double b1, double b2 )
{
	return sqrt(var/(2*pi)) * (std::exp(-std::pow(b1-mean,2)/2/var) - std::exp(-std::pow(b2-mean,2)/2/var));
}

/*
 * Computes \int_{b1}^{b2} dist(x) dx
 */
double I0( double mean, double var, double b1, double b2 )
{
	return erf( (b2-mean)/sqrt(2*var) ) / 2 - erf( (b1-mean)/sqrt(2*var) ) / 2;
}

/*
 * Computes the conditional mean inside a quantization interval
 */
double cmean( double mean, double var, double b1, double b2 )
{
	double	i0 = I0( mean, var, b1, b2 );
	double	i1 = I1( mean, var, b1, b2 );;
	double	rv;
	if( i0 != 0 )
		rv = i1/i0;
	else if( b1 == -inf && b2 == inf )
		rv = 0;
	else if( b1 == -inf )
		rv = b2;
	else if( b2 == inf )
		rv = b1;
	else
		rv = (b1+b2)/2;

	return rv;
}

//===========================================================================
// Rate distortion functions
//===========================================================================
double distortion_rate_function( double R, const function<double(double)>& phi )
{
	// Rate as a function of the dummy variable \lambda
	auto dist = [&]( double w, double l) {return min( phi(w), l );};
	auto Dist = [&]( double l ) {
		function<double(double)> fun = [&](double w){return dist(w,l);};
		return 1/(2*pi) * integral( fun, -pi, pi );
	};

	// Distortion as a function of the dummy variable \lambda
	auto rate = [&]( double w, double l) {return max( 1./2 * log2(phi(w)/l), 0. );};
	auto Rate = [&]( double l ) {
		function<double(double)> fun = [&](double w){return rate(w,l);};
		return 1/(2*pi) * integral( fun, -pi, pi );
	};

	// Guess value for \lambda
	int N     = 1000;
	rvec W    = linspace( -pi, pi, N );
	rvec ph   = entrywise( phi )(W);

	// Find \lambda
	auto fun    = [&] (double l) {	return Rate(l)-R; };
	double lmin = max(1e-12,min(ph));
	double lmax = max(ph);
	double l    = fzero( fun, lmin, lmax, 1000, 1e-12, 1e-6 );

	// Check
	if( abs(R-Rate(l))/R > 1e-4 )
		error("distortion_rate_function: Precision not attained");

	// Return distortion
	double D = Dist(l);
//	cout << "Dist = " << 10*log10(100*D) << endl;
//	cout << "Rate = " << 100*Rate(l) << endl;
	return D;
}

//===========================================================================
// Linear predictive quanitization (LPQ)
//===========================================================================
//---------------------------------------------------------------------------
// Scalar
//---------------------------------------------------------------------------
///*
// * Compute predictor
// */
//rarma lpc_predictor( const rseq &fx, int N )
//{
//	return 	lpc_predictor( mat<rseq>({{fx}}), N )(0,0);
//}
//
///*
// * Compute quantizer
// */
//quantizer lpc_quantizer( int L, const rseq &p, const rseq &fx, int K )
//{
//	// Design vector quantizer of one dimension
//	sequence<rmat> P  = ms2sm( mat<rseq>{{p}}  );
//	sequence<rmat> Fx = ms2sm( mat<rseq>{{fx}} );
//	vector_quantizer VQ = lpc_quantizer( L, P, Fx, K );
//
//	// Build scalar quantizer
//	int N = VQ.dict.size();
//	rvec dict( N );
//	for( int n = 0; n < N; n++ )
//		dict(n) = VQ.dict(n)(0);
//	return quantizer( dict );
//}
//
///*
// * Quantize a scalar signal given the PREDICTOR and QUANTIZER
// */
//rseq lpc( const rseq &x, rarma p, const quantizer &q )
//{
//	int T  = x.size();
//	int t1 = x.t1();
//	rseq xq( zeros(T), t1 );						// Quantized signal
//	double xp = 0;									// Predicted value
//	for( int t = x.t1(); t <= x.t2(); t++ ) {
//		double pe = x(t) - xp;
//		xq(t)     = q( pe ) + xp;
//	    xp        = p( xq(t) );
//	}
//
//	return xq;
//}

//---------------------------------------------------------------------------
// Vector
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Compute predictor
//---------------------------------------------------------------------------
control::arma<rmat,rvec> lpc_predictor( const sequence<rmat> &Fx_, int N )
{
	mat<rseq> Fx = sm2ms(Fx_);

	// Compute R_ and M_
	mat<rseq> R = Fx * adjoint(Fx);
	int L   = R.size1();
	rmat R_ = zeros(N*L,L);
	rmat M_ = zeros(N*L,N*L);
	for( int n = 0; n < N; n++ ) {
	    int idn1 = n*L;
	    int idn2 = idn1+L-1;
	    R_( range(idn1,idn2+1), range(0,L) ) = sample_get( R, -n-1 );
	    for( int m = 0; m < N; m++ ) {
	        int idm1 = m*L;
	        int idm2 = idm1+L-1;
	        M_( range(idn1,idn2+1), range(idm1,idm2+1) ) = sample_get( R, m-n ) ;
	    }
	}

	// Compute P_
	rmat P_ = lls( M_, R_ );

	// Build numerators
	mat<rvec> B(L,L);
	for( int i = 0; i < L; i++ )
		for( int j = 0; j < L; j++ )
			B(i,j) = zeros(N);

	// Fill numerators
	for( int n = 0; n < N; n++ ) {
	    int idn1 = n*L;
	    int idn2 = idn1+L-1;
	    rmat tmp = trans( P_( range(idn1,idn2+1), range(0,P_.size2()) ) );
		for( int i = 0; i < L; i++ )
			for( int j = 0; j < L; j++ )
				B(i,j)(n) = tmp(i,j);
	}

	// Build arma model
	mat<control::transfer_function> P(L,L);
	for( int i = 0; i < L; i++ )
		for( int j = 0; j < L; j++ )
			P(i,j).set_coeffs( B(i,j), rvec({1}) );

	return mtf2tfm(P);
}

//---------------------------------------------------------------------------
// Compute the quantizer using the open loop algorithm
//---------------------------------------------------------------------------
vector_quantizer lpc_quantizer_openloop( int L, const function<sequence<rvec>(sequence<rvec>)>& P,
		const sequence<rmat> &Fx, int K )
{
	// Training signal
	sequence<rvec> E(K,0);
	int J = Fx(Fx.t1()).size2();
	for( int k = 0; k < K; k++ )
		E(k) = randn(J);
	sequence<rvec> X  = (Fx * E).trunc(0,K-1);
	sequence<rvec> Xp = time_shift( P(X), 1 );
	sequence<rvec> Xe = X - Xp;

	// Design quantizer
	vector_quantizer VQ( Xe.buffer(), L );

	return VQ;
}

//---------------------------------------------------------------------------
// Compute the quantizer using the closed loop algorithm in:
// H. Khalil, K. Rose, and S. Regunathan, “The asymptotic closed-loop
// approach to predictive vector quantizer design with application in video
// coding,” IEEE Trans. Image Process., vol. 10, no. 1, pp. 15–23, Jan. 2001.
//---------------------------------------------------------------------------
vector_quantizer lpc_quantizer_closedloop( int L, const function<sequence<rvec>(sequence<rvec>)>& P,
		const sequence<rmat> &Fx, int K )
{
	// Training signal
	sequence<rvec> E(K,0);
	int J = Fx(Fx.t1()).size2();
	for( int k = 0; k < K; k++ )
		E(k) = randn(J);
	sequence<rvec> X = (Fx * E).trunc(0,K-1);

	// Main loop
	sequence<rvec> Xq(K,0);									// Quantized signal
	int I = Fx(Fx.t1()).size1();
	for( int t = Xq.t1(); t <= Xq.t2(); t++ )
		Xq(t) = zeros(I);
	vector_quantizer VQ;								// Quantizer
	vector_quantizer VQ0;
	double err = inf;
	double err0;
	do {
		sequence<rvec> Xp__= P(Xq);
		sequence<rvec> Xp_ = time_shift( Xp__, 1 );
		sequence<rvec> Xp  = Xp_.trunc(Xq.t1()+1, Xq.t2());	// Prediction
		sequence<rvec> Xe  = X - Xp;							// Prediction error
		VQ0 = VQ;
		VQ.train( Xe.buffer(), L );						// Quantizer
		Xq = VQ( Xe ) + Xp;								// Quantization

		// stop criterion
		err0 = err;
		err  = 0;
		for( int t = X.t1(); t <= X.t2(); t++ )
			err += norm(X(t)-Xq(t));

	} while( err < err0 );

	return VQ0;
}

//---------------------------------------------------------------------------
// Quantize a vector signal given the PREDICTOR and QUANTIZER
//---------------------------------------------------------------------------
sequence<rvec> lpc( const sequence<rvec> &x, const function<rvec(rvec)>& p, const vector_quantizer &q )
{
	int M  = x(x.t1()).size();
	int T  = x.size();
	int t1 = x.t1();
	sequence<rvec> xq( T, t1 );							// Quantized signal
	rvec xp = zeros(M);								// Predicted value
	for( int t = x.t1(); t <= x.t2(); t++ ) {
		rvec pe = x(t) - xp;
		xq(t)   = q( pe ) + xp;
		xp      = p( xq(t) );
	}

	return xq;
}

//---------------------------------------------------------------------------
// Quantize a vector signal given the PREDICTOR LENGTH and QUANTIZATION LEVELS
//---------------------------------------------------------------------------
vec<rseq> lpc( const vec<rseq> &x, const mat<rseq> &Fx, int P, int L, int K,
		const string &quantization )
{
	// compute the predictor
	sequence<rmat> Fx_                  = ms2sm(Fx);
	control::arma<rmat,rvec> p = lpc_predictor( Fx_, P );

	// Quantizer design
	vector_quantizer q;
	char mode = quantization.data()[0];
	if( mode == 'o' || mode == 'O' )
		q = lpc_quantizer_openloop( L, p, Fx_, K );
	else if( mode == 'c' || mode == 'C' )
		q = lpc_quantizer_closedloop( L, p, Fx_, K );
	else
		error( string("quantization mode ") + quantization + string(" not supported") );

	// Linear predictive quantization
	sequence<rvec> x_  = vs2sv( x );
	sequence<rvec> xq  = lpc( x_, p, q );
	vec<rseq> xq_ = sv2vs(xq);

	return xq_;
}

}
}
