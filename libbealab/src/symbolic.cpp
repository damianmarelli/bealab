#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOSYMBOLIC

#include <bealab/core/prelim.hpp>
#include <bealab/scilib/symbolic.hpp>

namespace bealab
{
namespace symbolic
{
// Module initialization
__attribute__((constructor))
static void init()
{
	system("rm GiNaC*.so 2> /dev/null");
}

// Non-commutative symbolic object ---------------------------------------------
void symbol_nc::print( const print_dflt &c, unsigned level ) const
{
	c.s << name;
}

void symbol_nc::print_latex( const GiNaC::print_latex &c, unsigned level ) const
{
	c.s << latex_name;
}

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT( symbol_nc, basic,
       print_func<print_dflt>( &symbol_nc::print).
       print_func<GiNaC::print_latex>( &symbol_nc::print_latex) );

/// @cond NEVER
int symbol_nc::compare_same_type( const basic& other ) const
{
	const symbol_nc& o = static_cast<const symbol_nc&>(other);
	int cmpval = name.compare(o.name);
	if( cmpval == 0 )
		return 0;
	else if( cmpval < 0 )
		return -1;
	else
		return 1;
}
/// @endcond

// Maxima symbols --------------------------------------------------------------
symbol plus("plus");
symbol minus("minus");
symbol inf("inf","\\infty");
symbol minf("minf","-\\infty");
symbol und("und");
symbol ind("ind");
symbol infinity("infinity");

// Numeric evaluation-----------------------------------------------------------
template<>
double numeval( const ex& x )
{
	return ex_to<GiNaC::numeric>( x.evalf() ).to_double();
}

template<>
complex numeval( const ex& x )
{
	numeric num = ex_to<GiNaC::numeric>( x.evalf() );
	double real = num.real().to_double();
	double imag = num.imag().to_double();
	return {real,imag};
}

bealab::function<double(double)> makefun( const ex& fun, const symbol& x )
{
	FUNCP_CUBA fp;
	ex fun1 = fun.subs( Pi==pi );
	compile_ex( lst(fun1), lst(x), fp );
	return [fp]( double x )
	{
		int ndim  = 1;
		int ncomp = 1;
		double res;
		fp( &ndim, &x, &ncomp, &res );
		return res;
	};
}

//bealab::function<double(double)> makefun( const ex& fun, const symbol& x )
//{
//	return [&fun,&x]( double x_ )
//	{
//		return numeval( fun.subs(x==x_) );
//	};
//}

// Create a matrix with non-commutative symbols --------------------------------
matrix symbolic_matrix_nc( int I, int J, const string& var )
{
	matrix A(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ ) {
			ostringstream screen_name, latex_name;
			screen_name << var << i << j;
			latex_name << var << "_{" << i << j << "}";
			A(i,j) = symbol_nc( screen_name.str(), latex_name.str() );
		}
	return A;
}

// Function ERF ----------------------------------------------------------------
namespace erf_ns
{
ex eval( const ex& x )
{
	if( x == 0 )
		return 0;
	else if( x == inf )
		return 1;
	else if( x == minf )
		return -1;
	else
		return erf(x).hold();
}

ex evalf( const ex& x )
{
	if( is_a<numeric>(x) )
		return bealab::erf( ex_to<numeric>(x).to_double() );
	else
		return erf(x).hold();
}

ex deriv( const ex& x, unsigned diff_param )
{
	return 2 / sqrt(Pi) * exp( -pow(x,2) );
}

unsigned serial = function::register_new(
		function_options( "erf", 1 ).
		eval_func( eval ).
		evalf_func( evalf ).
		derivative_func( deriv ) );
}

const function erf( const ex&e )
{
	return function( erf_ns::serial, e );
}

// Function TRANSPOSE ----------------------------------------------------------
namespace trans_ns
{
ex eval( const ex& e )
{
	if( is_a<matrix>(e) )
		return trans( ex_to<matrix>(e) );
	else if( is_a<function>(e)
			&& ex_to<function>(e).get_name() == string("trans") )
		return ex_to<function>(e).op(0);
	else
		return trans(e).hold();
}

unsigned serial = function::register_new(
		function_options( "trans", 1 ).
		eval_func( eval ) );
}

const function trans( const ex&e )
{
	return function( trans_ns::serial, e );
}

// Calculus --------------------------------------------------------------------
ex integrate( const ex& f, const ex& x )
{
	return maxima::functor("integrate")( f, x );
}

ex integrate( const ex& f, const ex& x, const ex& a, const ex& b )
{
	return maxima::functor("integrate")( f, x, a, b );
}

ex limit( const ex& f, const ex& x, const ex& val )
{
	return maxima::functor("limit")( f, x, val );
}

// Harmonic analysis transforms ------------------------------------------------
ex fourier_transform( const ex& fun, const realsymbol& t, const realsymbol& f )
{
	return integrate( fun * exp(-I*2*Pi*f*t), t, -inf, inf );
}

ex ifourier_transform( const ex& fun, const realsymbol& t, const realsymbol& f )
{
	return integrate( fun * exp(I*2*Pi*f*t), f, -inf, inf );
}

ex convolution( const ex& f1, const ex& f2, const ex& x )
{
	symbol tau;
	return integrate( f1.subs(x==tau) * f2.subs(x==x-tau), tau, -inf, inf );
}

}
}
#endif
