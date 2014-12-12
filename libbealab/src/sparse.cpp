#include <bealab/extensions/signal/sparse.hpp>

namespace bealab
{
namespace sparse
{

//------------------------------------------------------------------------------
// Class relaxation
//------------------------------------------------------------------------------
void relaxation::approximate()
{
	// Initialize the zth_indexes
	if( zth_indexes.size() == 0 )
		zth_indexes = vrange( 0, guess.size() );

	// Run the approximation
	rvec coeffs = optimize();

	// Build the result
	this->indexes.resize(0);
	this->coefficients.resize(0);
	double M = max( coeffs );
	int N    = coeffs.size();
	for( int n = 0; n < N; n++ ) {

		// If the index is in zth_indexes and its value is very small,
		// ignore it
		if( find( zth_indexes == n ).size() &&
			abs(coeffs(n)) < zero_threshold * M )
			continue;

		// Add the index and coefficient to the result
		this->indexes      = ivec{ this->indexes, {n} };
		this->coefficients = rvec{ this->coefficients, {coeffs(n)} };
	}
}

//------------------------------------------------------------------------------
// Class lp_relaxation
//------------------------------------------------------------------------------
rvec lp_relaxation::optimize_direct()
{
	// Optimization problem
	int N = guess.size();
	optimization::base* popt;
	if( this->tolerance == 0 )
		popt = new optimization::sqp(N);
	else
		popt = new optimization::mma(N);
	optimization::base& opt  = *popt;

	// Objective
	auto objfun = [this]( const rvec& x )
	{
		double val = norm( x, this->p );
		if( this->trace )
			cout << "lp_relaxation - norm       = " << val << endl;
		return val;
	};
	auto objfun_grad = [this]( const rvec& x )
	{
		int N = x.size();
		rvec g(N);
		for( int n = 0; n < N; n++ ) {
			if( x(n) == 0 ) {
				g(n) = 0;
			}
			else {
				double y = this->p * real( pow( abs(x(n)), this->p-1 ) );
				g(n)     = x(n)>0 ? y : -y;
			}
		}
		return g;
	};
	opt.set_objective( objfun, objfun_grad );

	// Constraint representing the sparse approximation problem
	auto erfun = [&]( const rvec& x )
	{
		double val = this->errfun(x) - this->tolerance;
		if( this->trace )
			cout << "lp_relaxation - constraint = " << val << endl;
		return val;
	};
	auto grad = [&]( const rvec& x )
	{
		return this->errfun_gradient(x);
	};
	if( this->tolerance )
		opt.add_inequality_constraint( erfun, grad );
	else
		opt.add_equality_constraint( erfun, grad );

	// Guess
	opt.guess = guess.size() == 0 ? zeros(N) : guess;

	// Run the optimization
	opt.trace = this->trace;
	rvec rv   = opt.optimize();

	// Return
	delete popt;
	return rv;
}

rvec lp_relaxation::optimize_recast()
{
	// Optimization problem
	int N = guess.size();
	optimization::base* popt;
	if( this->tolerance == 0 )
		popt = new optimization::sqp(2*N);
	else
		popt = new optimization::sqp(2*N);
	optimization::base& opt  = *popt;

	// Objective
	figure fig;
	timer t;
	opt.set_objective(
		[&](const rvec& z)
		{
			auto b     = z( range( N, z.size() ) );
			rvec bp    = real( pow( b, this->p ) );
			double val = sum(bp );
			if( this->trace )
				cout << "lp_relaxation - objective  = " << val << endl;
			if( this->plot && t.elapsed() > 1) {
				rvec x = z( range(0,N) );
				fig.clear().
					plot( x ).
					plot( rvec(b), blue, dotted ).
					plot( rvec(-b), blue, dotted );
				t.reset();
			}
			return val;
		},
		[&](const rvec& z)
		{
			rvec gx = zeros(N);
			auto b  = z( range(N,z.size()) );
			rvec gb = this->p * real( pow( b, this->p-1 ) );
			return rvec{gx,gb};
		}
	);

	// Constraints resulting from norm minimization
	opt.add_inequality_vconstraint(
		[&](const rvec& z)
		{
			auto x = z( range(0,N) );
			auto b = z( range(N,z.size()) );
			return x - b;
		},
		[&](const rvec&)
		{
			return rmat{{eye(N),-eye(N)}};
		}
	);
	opt.add_inequality_vconstraint(
		[&](const rvec& z)
		{
			auto x = z( range(0,N) );
			auto b = z( range(N,z.size()) );
			return -x - b;
		},
		[&](const rvec&)
		{
			return rmat{{-eye(N),-eye(N)}};
		}
	);

	// Constraint representing the sparse approximation problem
	auto erfun = [&](const rvec& z)
	{
		auto x = z( range(0,N) );
		double val = this->errfun(x) - this->tolerance;
		if( this->trace )
			cout << "lp_relaxation - constraint = " << val << endl;
		return val;
	};
	auto grad  = [&](const rvec& z)
	{
		auto x  = z( range(0,N) );
		rvec gx = this->errfun_gradient(x);
		rvec gb = zeros(N);
		return rvec{gx,gb};
	};
	if( this->tolerance )
		opt.add_inequality_constraint( erfun, grad );
	else
		opt.add_equality_constraint( erfun, grad );

	// Lower bounds
	opt.lower_bounds.resize(2*N);
	rvec lb_x                        = -inf * ones(N);
	rvec lb_b                        = 1e-8 * ones(N);
	opt.lower_bounds( range(0,N) )   = lb_x;
	opt.lower_bounds( range(N,2*N) ) = lb_b;

	// Guess
	opt.guess.resize(2*N);
	rvec x_guess              = guess.size() == 0 ? zeros(N) : guess;
	rvec b_guess              = max(abs(x_guess)) * ones(N);
	opt.guess( range(0,N) )   = x_guess;
	opt.guess( range(N,2*N) ) = b_guess;

	// Run the optimization
	opt.trace = this->trace;
	rvec rv = opt.optimize()( range(0,N) );

	// Return
	delete popt;
	return rv;
}

rvec lp_relaxation::optimize()
{
	rvec coeffs;
	if( direct_method )
		coeffs = optimize_direct();
	else
		coeffs = optimize_recast();
	return coeffs;
}

//------------------------------------------------------------------------------
// Class lp_relaxation_vector
//------------------------------------------------------------------------------
lp_relaxation_vector::lp_relaxation_vector()
{
	// Approximation error function
	errfun = [this]( const rvec& x )
	{
		return pow( norm(this->target-this->matrix * x), 2 );
	};

	// Gradient of the approximation error function
	errfun_gradient = [this]( const rvec& x )
	{
		return 2 * ( this->gramian*x - this->innerprods );
	};
}

void lp_relaxation_vector::approximate()
{
	gramian    = trans(matrix) * matrix;
	innerprods = trans(matrix) * target;
	lp_relaxation::approximate();
}

//------------------------------------------------------------------------------
// Class shrinking_gaussian
//------------------------------------------------------------------------------
rvec shrinking_gaussian::optimize( double width, const rvec& guess )
{
	// init
	double epsilon = 0;

//	// Set the objective function
//	auto objfun = [this,width]( const rvec& x )
//	{
//		int N = x.size();
//		double val = 0;
//		for( int n = 0; n < N; n++ ) {
//			double xn = x(n);
//			val += 1 - exp( -pow(xn/width,2) / 2 );
//		}
//
//		if( trace >= 2 )
//			cout << "hyder_mahata - objective  = " << val << endl;
//		return val;
//	};
//	auto objfun_grad = [this,width]( const rvec& x )
//	{
//		rvec g = zeros(x.size());
//		int N = x.size();
//		for( int n = 0; n < N; n++ ) {
//			int idx   = n;
//			double xn = x(idx);
//			g(idx)    = xn / pow(width,2) * exp( -pow(xn/width,2) / 2 );
//		}
//		return g;
//	};
	// Set the objective function
	auto objfun = [this,width]( const rvec& x )
	{
		int N = zth_indexes.size();
		double val = 0;
		for( int n = 0; n < N; n++ ) {
			int idx   = zth_indexes(n);
			double xn = x(idx);
			val += 1 - exp( -pow(xn/width,2) / 2 );
		}

		if( trace > 1 )
			cout << "shrinking_gaussian: objective  = " << val << endl;
		return val;
	};
	auto objfun_grad = [this,width]( const rvec& x )
	{
		rvec g = zeros(x.size());
		int N = zth_indexes.size();
		for( int n = 0; n < N; n++ ) {
			int idx   = zth_indexes(n);
			double xn = x(idx);
			g(idx)    = xn / pow(width,2) * exp( -pow(xn/width,2) / 2 );
		}
		return g;
	};

	// Set the constraints representing the approximation error tolerance
	auto constr = [this,&epsilon]( const rvec& x ) -> rvec
	{
		double val = errfun(x) - tolerance - epsilon;
		if( trace > 1 )
			cout << "shrinking_gaussian: constraint = " << val << endl;
		return {val};
	};
	auto constr_grad = [this]( const rvec& x ) -> rmat
	{
		rmat jac(1,x.size());
		jac.row(0) = errfun_gradient(x);
		return jac;
	};

	// Define the optimization problem
	int N = guess.size();
//	optimization::var2 sopt( N );
//	sopt.stop_fincrement_relative = 1e-9;
//	optimization::barrier opt( N, sopt );
//	opt.stop_duality_gap = 1e-9;
	optimization::mma opt( N );
	opt.stop_fincrement_relative = 1e-6;
	opt.stop_xincrement_relative = 1e-6;
	opt.set_objective( objfun, objfun_grad );
	opt.add_inequality_vconstraint( constr, constr_grad );
	opt.guess = guess;
	opt.trace = this->trace - 2;

	// Make the problem feasible
	epsilon = max( constr(guess) ) + 1e-6;
	epsilon = max( epsilon, 0. );
//	opt.stop_constraint_tolerance = 1e-8 + e;

	// Run the optimization
	rvec rv = opt.optimize();

	// Check if the soluiton is feasible
	bool c = sum( constr(rv) > opt.stop_constraint_tolerance );
	if( c ) {
//		rv = guess;
		if( trace )
			cout << "shrinking_gaussian: the optimization algorithm failed! ***" << endl;
	}

	return rv;
}

rvec shrinking_gaussian::optimize()
{
	// Initialization
	double width = max( abs(guess(indirect(zth_indexes))) );
	rvec x       = guess;

	// Main loop
	for(;; width *= .5 ) {

		// Update
		rvec x0     = x;
		x           = optimize( width, x0 );
		double rinc = norm(x-x0) / norm(x0);

		// Trace
		if( trace > 0 )
			cout << "shrinking_gaussian: width = " << width
				 << ", relative parameter increment = " << rinc << endl;

		// Stopping condition
		double M = max(abs(x(indirect(zth_indexes))));
		if( width/M <  zero_threshold/5 )
			break;
	}
	return x;
}

//------------------------------------------------------------------------------
// Class shrinking_gaussian
//------------------------------------------------------------------------------
//rvec shrinking_gaussian::optimize( double t, double gw, const rvec& wstart )
//{
//	// Choose the optimization algorithm
//	int N = wstart.size();
//	optimization::var2 opt( N );
//	opt.stop_fincrement_relative = 1e-9;
//
//	// Set the objective function
//	auto objfun = [this,gw,t]( const rvec& x )
//	{
//		// Gaussian weighting
//		int N = zth_indexes.size();
//		double f = 0;
//		for( int n = 0; n < N; n++ ) {
//			int idx   = zth_indexes(n);
//			double xn = x(idx);
//			f        += 1 - exp( -pow(xn/gw,2) / 2 );
//		}
//
//		// Barrier function
//		double ph = -log( tolerance - errfun(x) );
//
//		// Trace
//		if( this->trace ) {
//			cout << "shrinking_gaussian - Gaussian component = " << f << endl;
//			cout << "shrinking_gaussian - Barrier function   = " << ph << endl;
//		}
//
//		// Return the function value
//		return t*f + ph;
//	};
//	auto objfun_grad = [this,gw,t]( const rvec& x )
//	{
//		return t * grad_g(x,gw) + grad_b(x);
//	};
//	opt.set_objective( objfun, objfun_grad );
//
//	// Run the optimization
//	opt.guess = wstart;
//	opt.trace = this->trace;
//	return opt.optimize();
//}
//
//rvec shrinking_gaussian::optimize()
//{
//	// Initialization
//	rvec wstart = guess;
//	double gw   = max( abs(guess) );
//	rvec gg     = grad_g( guess, gw );
//	rvec gb     = grad_b( guess );
//	double t    = inner_prod( gb, gg ) / inner_prod( gg, gg );
//
//	// Main loop
//	for(;;) {
//
//		// Do one optimization
//		wstart = optimize( t, gw, wstart );
//
//		// Trace
//		if( trace )
//			cout << "shrinking_gaussian - gw = " << gw <<
//					", t = " << t << endl;
//
//		// Stopping condition
//		double M = max(abs(wstart));
//		if( gw/M < zero_threshold/5 )
//			break;
//
//		// Update
//		t  *= 2;
//		gw /= 2;
//	}
//
//	// Return the last warm-start vector
//	return wstart;
//}

}
}
