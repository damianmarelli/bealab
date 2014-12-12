#include <bealab/core/gsl.hpp>
#include <bealab/scilib/linalg.hpp>
#include <bealab/scilib/optimization.hpp>
#include <bealab/scilib/calculus.hpp>
#include <bealab/scilib/stats.hpp>
#include <gsl/gsl_multimin.h>
#include <nlopt.h>

namespace bealab
{
namespace optimization
{
//==============================================================================================
// Base class
//==============================================================================================
void base::set_objective( const function<double(const rvec&)>& fun,
					const function<rvec(const rvec&)>& grad )
{
	// Set the objective function
	objective.function = fun;

	// Choose the gradient
	function<rvec(const rvec&)> gradient;
	if( grad == 0 )
		gradient = bealab::gradient( fun );
	else
		gradient = grad;

	objective.gradient = gradient;
}

void base::add_inequality_constraint( const function<double(const rvec&)>& fun,
					const function<rvec(const rvec&)>& grad )
{
	inequality_constraints.push_back( fungrad() );
	inequality_constraints.back().function = fun;
	if( grad == 0 )
		inequality_constraints.back().gradient = bealab::gradient( fun );
	else
		inequality_constraints.back().gradient = grad;
}

void base::add_equality_constraint( const function<double(const rvec&)>& fun,
					const function<rvec(const rvec&)>& grad )
{
	equality_constraints.push_back( fungrad() );
	equality_constraints.back().function = fun;
	if( grad == 0 )
		equality_constraints.back().gradient = bealab::gradient( fun );
	else
		equality_constraints.back().gradient = grad;
}

void base::add_inequality_vconstraint( const function<rvec(const rvec&)>& fun,
					const function<rmat(const rvec&)>& jacob )
{
	inequality_vconstraints.push_back( vfunjac() );
	inequality_vconstraints.back().function = fun;
	if( jacob == 0 )
		inequality_vconstraints.back().jacobian = jacobian( fun );
	else
		inequality_vconstraints.back().jacobian = jacob;
}

void base::add_equality_vconstraint( const function<rvec(const rvec&)>& fun,
					const function<rmat(const rvec&)>& jacob )
{
	equality_vconstraints.push_back( vfunjac() );
	equality_vconstraints.back().function = fun;
	if( jacob == 0 )
		equality_vconstraints.back().jacobian = jacobian( fun );
	else
		equality_vconstraints.back().jacobian = jacob;
}

void base::add_inequality_vsconstraint( const function<rvec(const rvec&)>& fun,
					const function<rsmat(const rvec&)>& jacob )
{
	inequality_vsconstraints.push_back( vsfunjac() );
	inequality_vsconstraints.back().function = fun;
	inequality_vsconstraints.back().jacobian = jacob;
}

void base::add_equality_vsconstraint( const function<rvec(const rvec&)>& fun,
					const function<rsmat(const rvec&)>& jacob )
{
	equality_vsconstraints.push_back( vsfunjac() );
	equality_vsconstraints.back().function = fun;
	equality_vsconstraints.back().jacobian = jacob;
}

void base::test_gradient( const fungrad& fg, const rvec& x0 )
{
	rvec g  = fg.gradient(x0);
	rvec gn = bealab::gradient( fg.function, x0 );
	cout << "Provided gradient = " << g << endl;
	cout << "Numeric gradient  = " << gn << endl;
	cout << "Difference        = " << g - gn << endl;
	cout << "Ratio             = " << element_div( g, gn ) << endl;
	cout << "Maximum error     = " << max(abs(g-gn)) << endl;
	cout << "Relative error    = " << norm(g-gn) / norm(g) << endl;
}

void base::test_jacobian( const vfunjac& fj, const rvec& x0 )
{
	rmat J  = fj.jacobian(x0);
	rmat Jn = bealab::jacobian( fj.function, x0 );
	cout << "Provided Jacobian = " << endl;
	cout << J << endl;
	cout << "Numeric Jacobian  = " << endl;
	cout << Jn << endl;
	cout << "Difference        = " << endl;
	cout << J - Jn << endl;
	cout << "Ratio             = " << endl;
	cout << element_div( J, Jn ) << endl;
	cout << "Maximum error     = " << max(abs(J-Jn)) << endl;
	cout << "Relative error    = " << norm(J-Jn) / norm(J) << endl;
}

void base::test_derivatives( const rvec& x0 )
{
	// Objective
	cout << "Testing derivatives..." << endl;
	cout << "objective function" << endl;
	test_gradient( objective, x0 );

	// Scalar inequality constraints
	for( uint i = 0; i < inequality_constraints.size(); i++ ) {
		cout << "scalar inequality constraint #" << i << endl;
		test_gradient( inequality_constraints[i], x0 );
	}

	// Scalar equality constraints
	for( uint i = 0; i < equality_constraints.size(); i++ ) {
		cout << "scalar equality constraint #" << i << endl;
		test_gradient( equality_constraints[i], x0 );
	}

	// Vector inequality constraints
	for( uint i = 0; i < inequality_vconstraints.size(); i++ ) {
		cout << "Vector inequality constraint #" << i << endl;
		test_jacobian( inequality_vconstraints[i], x0 );
	}

	// Vector equality constraints
	for( uint i = 0; i < equality_vconstraints.size(); i++ ) {
		cout << "Vector equality constraint #" << i << endl;
		test_jacobian( equality_vconstraints[i], x0 );
	}
	cout << "... done testing derivatives..." << endl;
}

//==============================================================================================
// Quasi-Newton class
//==============================================================================================
//rmat quasinewton::inverse_hessian( const rmat& B, const rvec& x, const rvec& x0,
//		const rvec& g, const rvec& g0 )
//{
//	rvec s   = x - x0;
//	rvec q   = g - g0;
//	rvec v   = B * q;
//	rmat A   = outer_prod( v, s );
//	double c = inner_prod(s,q);
//	return B + ( 1 + inner_prod(q,v) / c ) * outer_prod(s,s) / c
//			 - ( trans(A) + A ) / c;
//}
rmat quasinewton::inverse_hessian( const rmat& B, const rvec& x, const rvec& x0,
		const rvec& g, const rvec& g0 )
{
	rvec s   = x - x0;
	rvec q   = g - g0;
	rvec v   = B * q;
	auto A   = ublas::outer_prod( v, s );
	double c = inner_prod(s,q);
	double d = sqrt( 1 + inner_prod(q,v) / c );
	rvec ds  = d * s;
	return B + ( ublas::outer_prod(ds,ds) - ( ublas::trans(A) + A ) ) / c;
}

rvec quasinewton::linesearch( const rvec& x0, const rvec& g, const rmat& B )
{
	// Init
	rvec x = x0;
	f0     = objective.function(x0);
	f      = f0;
	rvec d = - B * g;

	// Check if d is nan
	if( isnan(norm(d)) ) {
		if( trace >= 3 )
			cout << "backtracking: the search direction is nan ***" << endl;
		return x;
	}

	// Check if the search direction is descent. If not, reverse the search
	// direction.
	double drate = inner_prod(d,g);
	if( drate > 0 ) {
		if( trace >= 3 )
			cout << "backtracking: the descent rate in the search direction is positive = "
				 << drate << " ***" << endl;
		d = -d;
		drate = inner_prod(d,g);
	}

	// Trace
	if( trace >= 3 )
		cout << "backtracking: descent rate in the search direction = "  << drate << endl;

	// Main loop
	double t = 1;
	for(;; t *= beta ) {

		// Update
		rvec xd   = t * d;
		rvec xt   = x0 + xd;
		double ft = objective.function(xt);

		// Trace
		double dx = norm(xt-x0);
		if( trace >= 4 )
			cout << "backtracking: f(xt) = " << ft
				 << ", f(x0) = " << f0
				 << ", x-increment = " << dx << endl;

		// Stopping condition
		if( ft < f0 + alpha * t * drate || dx == 0 ) {
			x = xt;
			f = ft;
			break;
		}
	}

	return x;
}

bool quasinewton::stopping_condition( const rvec& x, const rvec& x0, const rvec& g )
{
	// Compute the conditions
	double xrinc = norm(x-x0) / norm(x0);
	double gn    = norm(g);
	double fabs  = f0-f;
	double frel  = fabs / abs(f0);

	// Trace
	if( trace >= 3 ) {
		cout << "quasinewton: relative increment in the parameters = " << xrinc << endl;
		cout << "quasinewton: gradient = " << gn << endl;
		cout << "quasinewton: relative decrement in the function = " << frel << endl;
		cout << "quasinewton: absolute decrement in the function = " << fabs << endl;
	}

	// Relative X-increment condition
	if( xrinc <= stop_xincrement_relative ) {
		if( trace >= 2 )
			cout << "quasinewton: stopping because of relative increment in the parameters = " << xrinc << endl;
		return true;
	}

	// Gradient condition
	if( gn <= stop_gradient ) {
		if( trace >= 2 )
			cout << "quasinewton: stopping because of gradient = " << gn << endl;
		return true;
	}

	// F-increment conditions
	if( frel <= stop_fincrement_relative ) {
		if( trace >= 2 )
			cout << "quasinewton: stopping because of relative decrement in the function = " << frel << endl;
		return true;
	}
	if( fabs <= stop_fincrement_absolute ) {
		if( trace >= 2 )
			cout << "quasinewton: stopping because of absolute decrement in the function = " << fabs << endl;
		return true;
	}

	return false;
}

rvec quasinewton::optimize()
{
	// Initialization
	int N  = guess.size();
	rvec x = guess;
//	while(true) {
//		if( objective.function(guess) == inf )
//			guess += eps* randn(guess.size());
//		else
//			break;
//	}
//	int count = 0;
//	while(true) {
//		if( objective.function(guess) == inf ) {
//			guess += eps * max(abs(guess)) * randn(guess.size());
//			if( count++ == 10 )
//				return guess;
//		}
//		else
//			break;
//	}
	rvec g = objective.gradient(guess);
	rmat B = eye(N);

	// Main loop
	while(true) {

		// Memorize previous values
		rvec x0   = x;
		rvec g0   = g;

		// Line search
		x = linesearch( x, g, B );
		g = objective.gradient(x);

		// Call the callback function
		if( callback_function )
			callback_function(x);

		if( trace >= 1 )
			cout << "quasinewton: f(x) = " << f << endl;

		// Update the inverse-Hessian
		B = inverse_hessian( B, x, x0, g, g0 );

//		//XXX
//		if( isnan(norm(B)) ) {
//			cout << "*** The inverse Hessian is nan ***" << endl;
//			cout << "Delta X    = " << norm(x-x0) << endl;
//			cout << "Delta grad = " << norm(g-g0) << endl;
//		}
//		else {
//			rvec e    = real( diag( get<0>( eig(B) ) ) );
//			double em = min( e );
//			if( em <= 0 ) {
//				cout << "*** The inverse Hessian has a non-positive eigenvalue. eigs = " << e << endl;
//				cout << "Delta X    = " << norm(x-x0) << endl;
//				cout << "Delta grad = " << norm(g-g0) << endl;
//			}
//		}

		// Stopping condition
		if( stopping_condition( x, x0, g ) )
			break;
	}
	return x;
}

//==============================================================================================
// Barrier class
//==============================================================================================
double barrier::barrier_function( const rvec& x )
{
	// Return value
	double val = 0;

	// Dense constraints
	int A = inequality_vconstraints.size();
	for( int i = 0; i < A; i++ ) {
		rvec c = inequality_vconstraints[i].function(x);
		c     -= stop_constraint_tolerance * ones(c.size());
		if( sum(c >= 0) )
			return inf;
		val += -sum( log( -c ) );
	}

	// Sparse constraints
	int B = inequality_vsconstraints.size();
	for( int i = 0; i < B; i++ ) {
		rvec c = inequality_vsconstraints[i].function(x);
		c     -= stop_constraint_tolerance * ones(c.size());
		if( sum(c >= 0) )
			return inf;
		val += -sum( log( -c ) );
	}
	return val;
}

rvec barrier::barrier_gradient( const rvec& x )
{
	// Return value
	int I = x.size();
	rvec g = zeros(I);

	// Dense constraints
	int A = inequality_vconstraints.size();
	for( int i = 0; i < A; i++ ) {
		rvec c = inequality_vconstraints[i].function(x);
		c     -= stop_constraint_tolerance * ones(c.size());
		int N = c.size();
		if( sum(c >= 0) )
			return -inf * ones(I);
		rmat J = inequality_vconstraints[i].jacobian(x);
		for( int n = 0; n < N; n++ )
			g -= J.row(n) / c(n);
	}

	// Sparse constraints
	int B = inequality_vsconstraints.size();
	for( int i = 0; i < B; i++ ) {
		rvec c = inequality_vsconstraints[i].function(x);
		c     -= stop_constraint_tolerance * ones(c.size());
		int N = c.size();
		if( sum(c >= 0) )
			return -inf * ones(I);
		rsmat J = inequality_vsconstraints[i].jacobian(x);
		for( int n = 0; n < N; n++ )
//			g -= J.row(n) / c(n);
			g -= ublas::matrix_row<rsmat>(J,n) / c(n);
	}

	return g;
}

rvec barrier::centering( const rvec& wstart )
{
	// Initialization of the centering algorithm
	cen.guess = wstart;

	// Objective of the centering algorithm
	auto objfun = [this]( const rvec& x )
	{
		double f  = objective.function(x);
		double ph = barrier_function(x);
		return t*f + ph;
	};
	auto objfun_grad = [this]( const rvec& x )
	{
		rvec grad_o = objective.gradient(x);
		rvec grad_b = barrier_gradient( x );
		return t * grad_o + grad_b;
	};
	cen.set_objective( objfun, objfun_grad );

	// Return the result of the centering algorithm
	return cen.optimize();
}

rvec barrier::optimize()
{
	// Checks
	assert( inequality_vconstraints.size() == 1 );

	// Initialization
	rvec wstart = guess;
	double M    = 0;
	int A = inequality_vconstraints.size();
	for( int i = 0; i < A; i++ )
		M += inequality_vconstraints[i].function(guess).size();
	int B = inequality_vsconstraints.size();
	for( int i = 0; i < B; i++ )
		M += inequality_vsconstraints[i].function(guess).size();
//	rvec go     = objective.gradient( guess );
//	rvec gb     = barrier_gradient( guess );
//	t           = abs( inner_prod( gb, go ) / inner_prod( go, go ) );
	t           = abs( barrier_function(guess) / objective.function(guess) );

	// Trace
	if( trace )
		cout << "Barrier IP method: f(x0) = "
			 << objective.function(wstart) << endl;

	// Main loop
	for(;; t *= 10 ) {

		// Solve the unconstrained sub-problem
		rvec x = centering( wstart );

		// See if to take the new value
		if( objective.function(x) <= objective.function(wstart) )
			wstart = x;

		// Call the callback function
		if( callback_function )
			callback_function(x);

		// Trace
		if( trace )
			cout << "Barrier IP method: t = " << std::setw(12) << t
				 << ", f(xt) = " << objective.function(wstart) << endl;

		// Stopping condition
		if( stop_duality_gap > M/t )
			break;
	}

	// Return the last warm-start vector
	return wstart;
}

//==============================================================================================
// GSL-based classes
//==============================================================================================
double gsl::objective_proxy( const void* x, void* par )
{
	return ((functors*)par)->fun( _gsl::vector( x ) );
}

void gsl::gradient_proxy( const void* x, void* par, void* gr )
{
	_gsl::vector g(gr);
	g = ((functors*)par)->grad( _gsl::vector( x ) );
}

void gsl::objgrad_proxy( const void* x, void* par, double* f, void* gr )
{
	rvec g;
	*f = ((functors*)par)->fungrad( _gsl::vector( x ), g );
	_gsl::vector gm( gr );
	gm = g;
}

bool gsl::check_stopping_condition()
{
	// Stop if error
	if( status ) {
		stopreason = error;
		if(trace) cerr << "GSLOPT: Stopping because of GSL error - " << gsl_strerror(status) << endl;
		if(trace) cerr << "\tgradient norm = " << norm(grad) << endl;
		return true;
	}

	// Stop if gradient
	if( norm(grad) <= stop_gradient ) {
		stopreason = gradient;
		if(trace) cout << "GSLOPT: Stopping because of gradient condition. Gradient norm = " << norm(grad) << endl;
		return true;
	}

	// Stop if fval
	if( fx < stop_fvalue ) {
		stopreason = fvalue;
		if(trace) cout << "GSLOPT: Stopping because of function value. f(x) = " << fx << endl;
		return true;
	}

	// Stop if absolute function increment
	if( (fx0-fx) < stop_fincrement_absolute ) {
		stopreason = fincrement_absolute;
		if(trace) cout << "GSLOPT: Stopping because of absolute function increment: f(x(t-1))-f(x(t)) = " << fx0-fx << endl;
		return true;
	}

	// Stop if relative function increment
	if( (fx0-fx)/fx0 < stop_fincrement_relative ) {
		stopreason = fincrement_relative;
		if(trace) cout << "GSLOPT: Stopping because of relative function increment: [f(x(t-1))-f(x(t))] / f(x(t-1)) = "
				       << (fx0-fx)/fx0 << endl;
		return true;
	}

//	// Stop if absolute parameter increment
//	if( norm(x0-x) < stop_xincrement_absolute ) {
//		stopreason = xincrement_absolute;
//		if(trace) cout << "GSLOPT: Stopping because of absolute parameter increment" << endl;
//		return true;
//	}

	// Stop if relative parameter increment
	if( norm(x0-x)/norm(x0) < stop_xincrement_relative ) {
		stopreason = xincrement_relative;
		if(trace) cout << "GSLOPT: Stopping because of relative parameter increment: ||x(t-1)-x(t)|| / ||x(t-1)|| = "
				       << norm(x0-x)/norm(x0) << endl;
		return true;
	}

	// Stop number of function evaluations
	if( Nfeval > stop_feval ) {
		stopreason = xincrement_relative;
		if(trace) cout << "GSLOPT: Stopping because of relative parameter increment: ||x(t-1)-x(t)|| / ||x(t-1)|| = "
				       << norm(x0-x)/norm(x0) << endl;
		return true;
	}

	return false;
}

rvec gsl_noder::optimize()
{
	// Callback parameters
	functors par = { objective.function, 0, 0 };

	// Prepare the objective function
	gsl_multimin_function objfunc;
	objfunc.n      = guess.size();
	objfunc.f      = reinterpret_cast<double(*)(const gsl_vector*,void*)>( objective_proxy );
	objfunc.params = &par;

	// Set minimization algorithm
	const gsl_multimin_fminimizer_type *T = reinterpret_cast<gsl_multimin_fminimizer_type*>(pproblem);
	gsl_multimin_fminimizer *s            = gsl_multimin_fminimizer_alloc( T, guess.size() );
	gsl_multimin_fminimizer_set( s, &objfunc, _gsl::vector(guess), _gsl::vector(.1*guess) );

	// Main loop
	size_t iter = 0;
	int status;
	do {
		// Do one iteration
		status = gsl_multimin_fminimizer_iterate( s );
		if( status ) {
			cerr<< gsl_strerror(status) << endl;
			break;
		}

		// Stop criterion
		double size = gsl_multimin_fminimizer_size( s );
		status      = gsl_multimin_test_size( size, 1e-6 );

		// Loop control
		iter++;

		// Call the callback function
		if( callback_function )
			callback_function(x);

	} while( status == GSL_CONTINUE && iter < 1e6 );

	// Finish
	rvec rv = _gsl::vector( s->x );
	gsl_multimin_fminimizer_free (s);
	return rv;
}

rvec gsl_der::optimize()
{
	// Functor for joint computation of function and gradient
	auto objgrad = [this](const rvec &x, rvec &g) {
		g = this->objective.gradient( x );
		return this->objective.function( x );
	};

	// Parameters with the functors to pass to the proxies
	functors par = { objective.function, objective.gradient, objgrad };

	// Prepare the objective function
	gsl_multimin_function_fdf objfunc;
	objfunc.n      = dim;
	objfunc.f      = reinterpret_cast<double(*)(const gsl_vector*,void*)>( objective_proxy );
	objfunc.df     = reinterpret_cast<void(*)(const gsl_vector*,void*,gsl_vector*)>( gradient_proxy );
	objfunc.fdf    = reinterpret_cast<void(*)(const gsl_vector*,void*,double*,gsl_vector*)>( objgrad_proxy );
	objfunc.params = &par;

	// Set minimization algorithm
	const gsl_multimin_fdfminimizer_type *algorithm = reinterpret_cast<gsl_multimin_fdfminimizer_type*>(pproblem);
	gsl_multimin_fdfminimizer *s                    = gsl_multimin_fdfminimizer_alloc( algorithm, dim );
	double tol       = 0.1;
	double step_size = 1;
	gsl_multimin_fdfminimizer_set( s, &objfunc, _gsl::vector(guess), step_size, tol );

	// Main loop
	fx     = inf;
	x      = inf * ones(dim);
	for( int i = 0;; i++ ) {
		fx0    = fx;
		x0     = x;
		status = gsl_multimin_fdfminimizer_iterate( s );
		x      = _gsl::vector( s->x );
		fx     = s->f;
		grad   = _gsl::vector( s->gradient );
		Nfeval = i;

		// Test gradient
		if( test_gradient ) {
			rvec ngrad = bealab::gradient(objective.function, x);
			cout << "Maximum gradient difference: " << max(abs(grad-ngrad)) << endl;
			if( trace ) {
				cout << "Provided gradient: "    << grad << endl;
				cout << "Numeric gradient:  "    << ngrad << endl;
			}
		}

		// Call the callback function
		if( callback_function )
			callback_function(x);

		// Trace
		if( trace )
			cout << "gsl_der: f(x) = " << fx << endl;

		// Stopping condition
		if( check_stopping_condition() )
			break;
	}

	// Return the optimum vector
	gsl_multimin_fdfminimizer_free( s );
	return x;
}

simplex::simplex( int dim ) : gsl_noder(dim)
{
	pproblem = const_cast<gsl_multimin_fminimizer_type*>(gsl_multimin_fminimizer_nmsimplex2);
}

conjugate_fr::conjugate_fr( int dim ) : gsl_der(dim)
{
	pproblem = const_cast<gsl_multimin_fdfminimizer_type*>(gsl_multimin_fdfminimizer_conjugate_fr);
}

conjugate_pr::conjugate_pr( int dim ) : gsl_der(dim)
{
	pproblem = const_cast<gsl_multimin_fdfminimizer_type*>(gsl_multimin_fdfminimizer_conjugate_pr);
}

bfgs::bfgs( int dim ) : gsl_der(dim)
{
	pproblem = const_cast<gsl_multimin_fdfminimizer_type*>(gsl_multimin_fdfminimizer_vector_bfgs2);
}

//==============================================================================================
// NLopt-based classes
//==============================================================================================
double nlopt::fungrad_proxy( unsigned N, const double* px, double* pgrad, void* pfunctors )
{
	// Vector argument
	rvec x(N);
	for( int n = 0; n < (int)N; n++ )
		x(n) = px[n];

	// Function value
	double val = reinterpret_cast<fungrad*>(pfunctors)->function( x );

	// Gradient
	if( pgrad ) {
		rvec grad = reinterpret_cast<fungrad*>(pfunctors)->gradient( x );
//		if( norm(grad) < stop_gradient )
//			nlopt_set_force_stop( *reinterpret_cast<nlopt_opt*>(pproblem), gradient );
		for( int n = 0; n < (int)N; n++ )
			pgrad[n] = grad(n);
	}

	// Return function value
	return val;
}

void nlopt::vfungrad_proxy( unsigned M, double* presult, unsigned N, const double* px, double* pgrad, void* pfunctors )
{
	// Vector argument
	rvec x(N);
	for( int n = 0; n < (int)N; n++ )
		x(n) = px[n];

	// Return function value
	rvec val = reinterpret_cast<vfunjac*>(pfunctors)->function( x );
	for( int m = 0; m < (int)M; m++ )
		presult[m] = val(m);

	// Gradient
	if( pgrad ) {
		rmat jacob = reinterpret_cast<vfunjac*>(pfunctors)->jacobian( x );
		for( int m = 0; m < (int)M; m++ )
			for( int n = 0; n < (int)N; n++ )
				pgrad[m*N+n] = jacob(m,n);
	}
}

nlopt::nlopt( int dim ) : base(dim)
{
	pproblem  = new nlopt_opt;
}

nlopt::~nlopt()
{
	if( pproblem ) {
		nlopt_destroy(*reinterpret_cast<nlopt_opt*>(pproblem));
		delete reinterpret_cast<nlopt_opt*>(pproblem);
	}
}

rvec nlopt::optimize()
{
	// Init
	nlopt_opt* pnlpprob = reinterpret_cast<nlopt_opt*>(pproblem);
	int C;
	double& tol = stop_constraint_tolerance;

	// Modify the objective function to include the callback function
	auto objfun = [this]( const rvec& x )
	{
		// Call the callback function
		if( callback_function )
			callback_function(x);

		// Return the objective function
		return objective.function(x);
	};

	// Set objective
	fungrad obj = { objfun, objective.gradient };
	nlopt_set_min_objective( *pnlpprob, fungrad_proxy, &obj );

	// Set scalar inequality constraints
	C = inequality_constraints.size();
	for( int c = 0; c < C; c++ )
		nlopt_add_inequality_constraint( *pnlpprob, fungrad_proxy,
				&inequality_constraints[c], tol );

	// Set scalar equality constraints
	C = equality_constraints.size();
	for( int c = 0; c < C; c++ )
		nlopt_add_equality_constraint( *pnlpprob, fungrad_proxy,
				&equality_constraints[c], tol );

	// Set vector inequality constraints
	C   = inequality_vconstraints.size();
	for( int c = 0; c < C; c++ ) {
		int nconstraints = inequality_vconstraints[c].function(zeros(dim)).size();
		rvec tol         = 1e-8 * ones(nconstraints);
		nlopt_add_inequality_mconstraint( *pnlpprob, nconstraints, vfungrad_proxy,
				&inequality_vconstraints[c], &tol(0) );
	}

	// Set vector equality constraints
	C   = equality_vconstraints.size();
	for( int c = 0; c < C; c++ ) {
		int nconstraints = inequality_vconstraints[c].function(zeros(dim)).size();
		rvec tol         = 1e-8 * ones(nconstraints);
		nlopt_add_equality_mconstraint( *pnlpprob, nconstraints, vfungrad_proxy,
				&equality_vconstraints[c], &tol(0) );
	}

	// Set bounds
	nlopt_set_lower_bounds( *pnlpprob, &lower_bounds(0) );
	nlopt_set_upper_bounds( *pnlpprob, &upper_bounds(0) );

	// Set stopping conditions
	nlopt_set_stopval( *pnlpprob, stop_fvalue );
	nlopt_set_ftol_rel( *pnlpprob, stop_fincrement_relative );
	nlopt_set_ftol_abs( *pnlpprob, stop_fincrement_absolute );
	nlopt_set_xtol_rel( *pnlpprob, stop_xincrement_relative );
	nlopt_set_xtol_abs( *pnlpprob, &stop_xincrement_absolute(0) );
	nlopt_set_maxeval( *pnlpprob, stop_feval );
	nlopt_set_maxtime( *pnlpprob, stop_time );

	// Optimize
	solution = guess;
	nlopt_result rv = nlopt_optimize( *pnlpprob, &solution(0), &minimum );

	// Set stop reason
	string reason;
	switch( rv ) {
	case NLOPT_SUCCESS:
		stopreason = none;
		reason     = "success";
		break;
	case NLOPT_STOPVAL_REACHED:
		stopreason = fvalue;
		reason     = "Stopping because of function value";
		break;
	case NLOPT_FTOL_REACHED:
		stopreason = fincrement;
		reason     = "Stopping because of function increment (relative or absolute)";
		break;
	case NLOPT_XTOL_REACHED:
		stopreason = xincrement;
		reason     = "Stopping because of parameter increment (relative or absolute)";
		break;
	case NLOPT_MAXEVAL_REACHED:
		stopreason = feval;
		reason     = "Stopping because of number of function evaluations";
		break;
	case NLOPT_MAXTIME_REACHED:
		stopreason = time;
		reason     = "Stopping because of time limit";
		break;
	case NLOPT_FAILURE:
		stopreason = error;
		reason     = "Stopping because of general failure";
		break;
	case NLOPT_INVALID_ARGS:
		stopreason = error;
		reason     = "invalid arguments";
		break;
	case NLOPT_OUT_OF_MEMORY:
		stopreason = error;
		reason     = "Stopping because of out of memory";
		break;
	case NLOPT_ROUNDOFF_LIMITED:
		stopreason = error;
		reason     = "Stopping because of roundoff error";
		break;
	case NLOPT_FORCED_STOP:
		stopreason = (stopcode)nlopt_get_force_stop( *pnlpprob );
		switch(stopreason) {
		case gradient:
			reason = "Stopping because of gradient";
			break;
		default:
			reason = "Stopping because of forced termination";
			break;
		}
		break;
	}

	if(trace) cout << "NLOPT: " << reason << endl;

	return solution;
}

neldermead::neldermead( int dim ) : nlopt(dim)
{
	*reinterpret_cast<nlopt_opt*>(pproblem) = nlopt_create( NLOPT_LN_NELDERMEAD, dim );
}

lbfgs::lbfgs( int dim ) : nlopt(dim)
{
	*reinterpret_cast<nlopt_opt*>(pproblem) = nlopt_create( NLOPT_LD_LBFGS, dim );
}

tnewton::tnewton( int dim ) : nlopt(dim)
{
	*reinterpret_cast<nlopt_opt*>(pproblem) = nlopt_create( NLOPT_LD_TNEWTON, dim );
}

tnewton_precond::tnewton_precond( int dim ) : nlopt(dim)
{
	*reinterpret_cast<nlopt_opt*>(pproblem) = nlopt_create( NLOPT_LD_TNEWTON_PRECOND, dim );
}

tnewton_restart::tnewton_restart( int dim ) : nlopt(dim)
{
	*reinterpret_cast<nlopt_opt*>(pproblem) = nlopt_create( NLOPT_LD_TNEWTON_RESTART, dim );
}

tnewton_precond_restart::tnewton_precond_restart( int dim ) : nlopt(dim)
{
	*reinterpret_cast<nlopt_opt*>(pproblem) = nlopt_create( NLOPT_LD_TNEWTON_PRECOND_RESTART, dim );
}

var1::var1( int dim ) : nlopt(dim)
{
	*reinterpret_cast<nlopt_opt*>(pproblem) = nlopt_create( NLOPT_LD_VAR1, dim );
}

var2::var2( int dim ) : nlopt(dim)
{
	*reinterpret_cast<nlopt_opt*>(pproblem) = nlopt_create( NLOPT_LD_VAR2, dim );
}

mma::mma( int dim ) : nlopt(dim)
{
	*reinterpret_cast<nlopt_opt*>(pproblem) = nlopt_create( NLOPT_LD_MMA, dim );
}

sqp::sqp( int dim ) : nlopt(dim)
{
	*reinterpret_cast<nlopt_opt*>(pproblem) = nlopt_create( NLOPT_LD_SLSQP, dim );
}

auglag::auglag( int dim, const nlopt& opt ) : nlopt(dim)
{
	nlopt_opt *popt = reinterpret_cast<nlopt_opt*>(pproblem);
	nlopt_opt *psub = reinterpret_cast<nlopt_opt*>(opt.get_pproblem());
	*popt = nlopt_create( NLOPT_LD_SLSQP, dim );
	nlopt_set_local_optimizer( *popt, *psub );
}

//==============================================================================================
// Ipopt-base primal-dual class
//==============================================================================================
#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOIPOPT

#include <coin/IpStdCInterface.h>

rvec primal_dual::ptr2vec( int N, double* px )
{
	rvec x(N);
	for( int n = 0; n < N; n++ )
		x(n) = px[n];
	return x;
}

void primal_dual::vec2ptr( const rvec& x, double* px )
{
	int N = x.size();
	for( int n = 0; n < N; n++ )
		px[n] = x(n);
}

void primal_dual::mat2ptr( const rmat& x, double* px )
{
	int I = x.size1();
	int J = x.size2();
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			px[i*J+j] = x(i,j);
}

int primal_dual::eval_f( int n, double* px, int new_x, double* obj_value,
		void* user_data )
{
	*obj_value = reinterpret_cast<primal_dual*>(user_data)->
			objective.function( ptr2vec(n,px) );
	return true;
}

int primal_dual::eval_grad_f( int n, double* px, int new_x, double* grad_f,
		void* user_data )
{
	vec2ptr( reinterpret_cast<primal_dual*>(user_data)->
			objective.gradient( ptr2vec(n,px) ),
			grad_f );
	return true;
}

int primal_dual::eval_g( int n, double* px, int new_x, int m, double* g,
		void* user_data )
{
	rvec grieq = reinterpret_cast<primal_dual*>(user_data)->
			inequality_vconstraints[0].function( ptr2vec(n,px) );
//	rvec greq  = reinterpret_cast<ip*>(user_data)->
//			equality_vconstraints[0].function( ptr2vec(n,px) );
//	rvec gr = { grieq, greq };
	rvec gr = grieq;
	vec2ptr( gr, g );
	return true;
}

int primal_dual::eval_jac_g( int J, double* px, int new_x,
		int I, int nele_jac, int* iRow, int *jCol,
		double* values, void* user_data )
{
	if( values == NULL )
		for( int i = 0; i < I; i++ )
			for( int j = 0; j < J; j++ ) {
				iRow[i*J+j] = i;
				jCol[i*J+j] = j;
			}
	else {
		rmat jacieq = reinterpret_cast<primal_dual*>(user_data)->
				inequality_vconstraints[0].jacobian( ptr2vec(J,px) );
//		rmat jaceq = reinterpret_cast<ip*>(user_data)->
//				equality_vconstraints[0].jacobian( ptr2vec(J,px) );
//		rmat jac = {{jacieq},{jaceq}};
//		mat2ptr( jac, values );
		mat2ptr( jacieq, values );
	}
	return true;
}

rvec primal_dual::optimize()
{
	// Take the number of constraints
	int Nieq = 0;
	if( inequality_vconstraints.size() )
		Nieq = inequality_vconstraints[0].function( guess ).size();
	int Neq = 0;
	if( equality_vconstraints.size() )
		Neq = equality_vconstraints[0].function( guess ).size();
	int N = Nieq + Neq;

	// Constraint bounds
	rvec lcb = { -1e19 * ones(Nieq), zeros(Neq) };
	rvec ucb = zeros(N);

	// Define the Ipopt problem
	IpoptProblem nlp = CreateIpoptProblem( dim, &lower_bounds(0),
			&upper_bounds(0), N, &lcb(0), &ucb(0), N*dim, 0, 0,
			eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h );

	// Set options
	AddIpoptNumOption( nlp, (char*)"tol", stop_fincrement_relative );
	AddIpoptIntOption( nlp, (char*)"print_level", 0 );
	AddIpoptStrOption( nlp, (char*)"mu_strategy", (char*)"adaptive" );
	AddIpoptStrOption( nlp, (char*)"hessian_approximation", (char*)"limited-memory" );

	// Solve the problem
	rvec x = guess;
//	ApplicationReturnStatus status =
	IpoptSolve( nlp, &x(0), NULL, NULL, NULL, NULL, NULL, this );

	// Free the allocated memory
	FreeIpoptProblem(nlp);

	return x;
}
#endif

}
}
