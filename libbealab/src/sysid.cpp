#include <bealab/scilib/fourier.hpp>
#include <bealab/scilib/optimization.hpp>
#include <bealab/extensions/control/sysid.hpp>

namespace bealab
{
namespace sysid
{
//------------------------------------------------------------------------------
// Class kalman_td
//------------------------------------------------------------------------------
void kalman_td::identify()
{
	// init
	const rvec &u = in;
	const rvec &y = out;
	int n         = max( n_num, n_den );
	int T         = u.size();
	int Tr        = T - n + 1;

	// Regresor matrix Phi
	rmat Phi_u( Tr, n_num );
	rmat Phi_y( Tr, n_den-1 );
	for( int tr = 0; tr < Tr; tr++ ) {
		int t            = tr + n - 1;
		Phi_u.row(tr) = flip( u( range(t-n_num+1,t+1) ) );
		Phi_y.row(tr) = -flip( y( range(t-n_den+1,t) ) );
	}
	rmat Phi = mconcat<double>{{Phi_u, Phi_y }};

	// Target vector yr
	rvec yr = y( range(n-1,y.size()) );

	// Solution (th)
	rvec th = lls( Phi, yr );
	num     = th( range(0,n_num) );
	den     = vconcat<double>{ {1}, th( range(n_num,th.size()) ) };
}

//------------------------------------------------------------------------------
// Class steiglitz_mcbride_td
//------------------------------------------------------------------------------
void steiglitz_mcbride_td::identify()
{
	// init
	const rvec &u = in;
	const rvec &y = out;

	// Main loop
	double	err = inf;
	double	derr;
	den = rvec{1};
	rvec num0, den0;
	do {

		// Previous estimation
		num0 = num;
		den0 = den;

		// Pre-filtering
		rvec uf = transfer_function( rvec{1}, den )( u ).vec();
		rvec yf = transfer_function( rvec{1}, den )( y ).vec();

		// Kalman ID
		kalman_td id( n_num, n_den );
		id.input(uf);
		id.output(yf);
		id.identify();
		num = id.num;
		den = id.den;

		// Stopping condition
		double erro = err;
		rvec yh     = transfer_function( num, den )( u ).vec();
		err         = norm(y-yh) / norm(y);
		derr        = err - erro;

	} while(derr < 0);

	// Keep the previous estimates
	num = num0;
	den = den0;
}

//------------------------------------------------------------------------------
// Class kalman_fd
//------------------------------------------------------------------------------
void kalman_fd::identify()
{
	// Compute the regressor matrix R
	int M = w.size();
	int N = n_num + n_den - 1;
	cmat R(M,N);
	for( int m = 0; m < M; m++ ) {
		for( int n = 0; n < n_num; n++ )
			R(m,n) = exp(-j*w(m)*n) * in(m);
		for( int n = 1; n < n_den; n++ )
			R(m,n+n_num-1) = - exp(-j*w(m)*n) * out(m);
	}

	// Convert to a real problem
	rmat Rr = {{real(R)},{imag(R)}};
	rvec yr = {{real(out)},{imag(out)}};

	// Estimate the parameters
	rvec par = lls( Rr, yr );
	num      = par(range(0,n_num));
	den      = { {1}, par(range(n_num,N))};
}

//------------------------------------------------------------------------------
// Class iterative_reweighting
//------------------------------------------------------------------------------
cvec iterative_reweighting::apply_fun( const cvec& x, const cvec& y )
{
	int N = w.size();
	cvec p(N);
	for( int n = 0; n < N; n++ )
		p(n)  = fun(x(n),y(n));
	return p;
}

void iterative_reweighting::identify()
{
	// Initialization
	kalman_fd id( n_num, n_den );
	id.input(in);
	id.output(out);
	id.frequencies(w);
	id.identify();
	num = id.num;
	den = id.den;

	// Main loop
	double	err = inf;
	rvec num0, den0;
	while(true) {

		// Pre-filtering
		cvec n  = dtft( rseq(num), w );											// Frequency response of the numerator
		cvec d  = dtft( rseq(den), w );											// Frequency response of the denominator
		cvec x  = element_prod( element_div( n, d ), in );						// n/d * in
		cvec p  = apply_fun( out, x );											// fun( out, n/d*in )
		cvec q  = element_prod(d,out) - element_prod(n,in);						// d*out - n*in
		cvec m  = element_div(p,q);												// ( log(out) - log(n/d*in) ) / ( d*out - n*in )
		cvec uf = element_prod( in, m );
		cvec yf = element_prod( out, m );

		// Stopping condition
		double erro = err;
		err         = pow( norm(p), 2 );
		if(trace) cout << "iterative_reweighting: Error = " << err << endl;
		if( err - erro >= 0 )
			break;

		// Remember previous estimation
		num0 = num;
		den0 = den;

		// Kalman ID
		kalman_fd id( n_num, n_den );
		id.input(uf);
		id.output(yf);
		id.frequencies(w);
		id.identify();
		num = id.num;
		den = id.den;
	}

	// Keep the previous estimates
	num = num0;
	den = den0;
}

//------------------------------------------------------------------------------
// Class steiglitz_mcbride_fd
//------------------------------------------------------------------------------
steiglitz_mcbride_fd::steiglitz_mcbride_fd( int nn, int dd ) :
		iterative_reweighting(nn,dd)
{
	fun = []( complex x, complex y ){ return x - y; };
}

//------------------------------------------------------------------------------
// Class kobayashi
//------------------------------------------------------------------------------
kobayashi::kobayashi( int nn, int dd ) :
		iterative_reweighting(nn,dd)
{
	fun = []( complex x, complex y ){ return log(x/y); };
}

//------------------------------------------------------------------------------
// Class quasinewton
//------------------------------------------------------------------------------
double quasinewton::objfun( const rvec& par )
{
	// Init
	int N    = n_num + n_den - 1;
	rvec num = par(range(0,n_num));
	rvec den = { {1}, par(range(n_num,N)) };
	cvec n   = dtft( rseq(num), w );
	cvec d   = dtft( rseq(den), w );
	cvec x   = element_prod( element_div( n, d ), in );

	// Sum of squares
	int I = w.size();
	rvec p(I);
	for( int i = 0; i < I; i++ )
		p(i) = pow( abs(fun(out(i),x(i))), 2 );
	return sum(p);
}

rvec quasinewton::objgrad( const rvec& par )
{
	// Init
	int N    = n_num + n_den - 1;
	rvec num = par(range(0,n_num));
	rvec den = { {1}, par(range(n_num,N)) };
	cvec n   = dtft( rseq(num), w );
	cvec d   = dtft( rseq(den), w );
	cvec x   = element_prod( element_div( n, d ), in );

	// Auxiliary weighting for the numerator and denominator
	int I = w.size();
	cvec an(I), ad(I);
	for( int i = 0; i < I; i++ ) {
		complex fdf = conj(fun(out(i),x(i))) * dfun(out(i),x(i));
		an(i)       = fdf * in(i) / d(i);
		ad(i)       = fdf * in(i) * n(i) / pow(d(i),2);
	}

	// Compute the gradient
	rvec g(N);
	for( int n = 0; n < n_num; n++ )
		g(n) =  2 * sum( element_prod( real(an), cos(-n*w) ) -
					     element_prod( imag(an), sin(-n*w) ) );
	for( int n = n_num; n < N; n++ )
		g(n) = -2 * sum( element_prod( real(ad), cos(-(n-n_num+1)*w) ) -
						 element_prod( imag(ad), sin(-(n-n_num+1)*w) ) );

	return g;
}

void quasinewton::identify()
{
	// Initialization
	if(trace) cout << "quasinewton_fd: initialization stage..." << endl;
	iterative_reweighting id( n_num, n_den );
	id.input(in);
	id.output(out);
	id.frequencies(w);
	id.fun   = fun;
	id.trace = max( trace-1, 0 );
	id.identify();

	// Optimization
	if(trace) cout << "quasinewton_fd: optimization stage..." << endl;
	int N = n_num + n_den - 1;
	optimization::quasinewton qn(N);
	qn.set_objective( [this](const rvec& x){ return objfun(x); },
					  [this](const rvec& x){ return objgrad(x); } );
	qn.guess = { id.num, id.den(range(1,n_den)) };
	qn.trace = max( trace-1, 0 );
//	qn.test_derivatives( qn.guess );
//	cin.get();
	rvec x   = qn.optimize();
	num      = x(range(0,n_num));
	den      = { {1}, x(range(n_num,N)) };
}

//------------------------------------------------------------------------------
// Class linear_fd
//------------------------------------------------------------------------------
linear_fd::linear_fd( int nn, int dd ) :
		quasinewton(nn,dd)
{
	fun   = []( complex x, complex y ){ return x-y; };
	dfun  = []( complex x, complex y ){ return -1; };
}

//------------------------------------------------------------------------------
// Class logarithmic_fd
//------------------------------------------------------------------------------
logarithmic_fd::logarithmic_fd( int nn, int dd ) :
		quasinewton(nn,dd)
{
	fun   = []( complex x, complex y ){ return log(x/y); };
	dfun  = []( complex x, complex y ){ return -1/y; };
}

/*
//------------------------------------------------------------------------------
// Class plinreg
//------------------------------------------------------------------------------
void plinreg::iteration()
{
	// Init
	const rvec &u = *p_in;
	const rvec &y = *p_out;
	int n         = max( n_num, n_den );
	int T         = u.size();
	int Tr        = T - n + 1;

	// Generate x
	rvec x = transfer_function( num, den )( u ).vec();

	// Regressor matrix (Phi)
	rmat Phi_u(Tr,n_num);
	rmat Phi_x(Tr,n_den-1);
	for( int tr = 0; tr < Tr; tr++ ) {
		int t = tr + n - 1;
		Phi_u.row( tr ) = flip( u( range(t-n_num+1,t+1) ) );
		Phi_x.row( tr ) = -flip( x( range(t-n_den+1,t) ) );
	}
	rmat Phi = {{ Phi_u, Phi_x }};

	// Target vector (yr)
	rvec yr = y( range(n-1,y.size()) );

	// th
	rvec th = lls( Phi, yr );
	num     = th( range(0,n_num) );
	den     = { {1}, th( range(n_num,th.size()) ) };
}

void plinreg::identify()
{
	// Init
	const rvec &u = *p_in;
	const rvec &y = *p_out;

	// Main loop
	double err = inf;
	double derr;
	rvec num0, den0;
	do {
		// Previous estimation
		num0 = num;
		den0 = den;

		// Run one iteration
		iteration();

		// Stopping condition
		rvec yh     = transfer_function( num, den )( u ).vec();
		double erro = err;
		err         = norm(y-yh) / norm(y);
		derr        = err - erro;

	} while(derr < 0);

	// Keep the previous estimate
	num = num0;
	den = den0;
}
*/

}
}
