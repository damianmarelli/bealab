/// @file bealab/core/python.hpp
/// Interface to call Python functions.

#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOPYTHON

#ifndef _BEALAB_PYTHON_
#define	_BEALAB_PYTHON_

#include <bealab/core/blas.hpp>

namespace bealab
{
namespace python
{

/// @defgroup interfaces_python Python interface
/// Interface to call Python functions.
/// @{

/// Interpreter class.
class object {
	void *pdata;
public:
	object() : pdata(NULL) {}
	object( void *p ) { pdata = p; }
	object( const object &o );
	object( int val );
	object( double val );
	object( complex val );
	object( const char *val );
	object( const string &val );
	template<class T>
	object( const Vec<T> &val );
	template<class T>
	object( const Mat<T> &val );
	~object();
	void operator=( const object &o );
	void* data() const { return pdata; }
	operator int() const;
	operator double() const;
	operator complex() const;
	operator string() const;
	template<class T>
	operator Vec<T>() const;
	template<class T>
	operator Mat<T>() const;
};

/// @name Internal functions
void _init();		///< Initializes the interface.
/// @}

/// Functor wrapper for python functions
class functor {
	void* pFunc;
	object _fcall( const deque<object> &O );
public:
	functor( const string& module, const string& function );
	template<class... T>
	object operator()( const T&... argin )
	{
	    // Parse parameters
	    deque<object> O;
	    O = parse_variadic_template<object>( argin... );

	    return _fcall( O );
	}
};

/// @name Plotting interface.
//static functor show( "pylab", "show" );
//static functor figure( "pylab", "figure" );
//static functor close( "pylab", "close" );
//static functor plot( "pylab", "plot" );
//static functor hold( "pylab", "hold" );
//static functor title( "pylab", "title" );
//static functor axes( "pylab", "axes" );
//static functor xlabel( "pylab", "xlabel" );
//static functor ylabel( "pylab", "ylabel" );
//static functor xlim( "pylab", "xlim" );
//static functor ylim( "pylab", "ylim" );
//static functor subplot( "pylab", "subplot" );
//static functor legend( "pylab", "legend" );
/// @}

/// @}
}
}
#endif

#endif
