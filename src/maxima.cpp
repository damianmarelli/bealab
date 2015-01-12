#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOSYMBOLIC

#include <bealab/core/prelim.hpp>
#include <bealab/core/maxima.hpp>
#include <ginac/ginac.h>

namespace bealab
{
namespace symbolic
{
namespace maxima
{

// Initialize the elements maxima::functor
const symbol functor::EULER__ = symbol("EULER__");
bool functor::trace = false;

// Function INTEGRATE2 ---------------------------------------------------------
namespace integrate2_ns
{
void print_latex_( const ex& fun, const ex& var, const print_context& c )
{
	c.s << "\\int " << fun << " d" << var;
}

unsigned serial = function::register_new(
			function_options( "integrate", 2 ).
			print_func<print_latex>( print_latex_ ).
			overloaded(2) );
}

//const function integrate2( const ex&e1, const ex&e2 )
//{
//	return function( integrate2_ns::serial, e1, e2 );
//}

// Function INTEGRATE4 ---------------------------------------------------------
namespace integrate4_ns
{
void print_latex_( const ex& fun, const ex& var, const ex& a,
		const ex& b, const print_context& c )
{
	c.s << "\\int_{" << a << "}^{" << b << "} " << fun << " d" << var;
}

unsigned serial = function::register_new(
			function_options( "integrate", 4 ).
			print_func<print_latex>(print_latex_).
			overloaded(2) );
}

//const function integrate4( const ex&e1, const ex&e2, const ex&e3, const ex&e4 )
//{
//	return function( integrate4_ns::serial, e1, e2, e3, e4 );
//}

// Function LIMIT --------------------------------------------------------------
namespace limit_ns
{
void print_latex_( const ex& f, const ex& x, const ex& val, const print_context& c )
{
	c.s << "\\lim_{" << x << "\\rightarrow" << val << "} " << f;
}

unsigned serial = function::register_new(
			function_options( "limit", 3 ).
//			eval_func( eval ).
			print_func<print_latex>(print_latex_) );
}

//const function limit( const ex& f, const ex& x, const ex& val )
//{
//	return function( limit_ns::serial, f, x, val );
//}

// Function DIFF3 --------------------------------------------------------------
namespace diff3_ns
{
ex eval( const ex& f, const ex& x, const ex& n )
{
	return diff( f, ex_to<symbol>(x), ex_to<numeric>(n).to_double() );
}

unsigned serial = function::register_new(
			function_options( "diff", 3 ).
			eval_func( eval ).
			overloaded(2) );
}

// Function DIFF5 --------------------------------------------------------------
namespace diff5_ns
{
ex eval( const ex& f, const ex& x, const ex& nx, const ex& y, const ex& ny )
{
	ex d1 = diff( f, ex_to<symbol>(x), ex_to<numeric>(nx).to_double() );
	return diff( d1, ex_to<symbol>(y), ex_to<numeric>(ny).to_double() );
}

unsigned serial = function::register_new(
			function_options( "diff", 5 ).
			eval_func( eval ).
			overloaded(2) );
}

}
}
}

#endif
