/// @file bealab/core/python.hpp
/// Interface to call Python functions.

#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOPYTHON

#ifndef _BEALAB_PYTHON_
#define	_BEALAB_PYTHON_

#include <bealab/core/blas.hpp>

namespace bealab
{
/// Python interface
namespace python
{

/// @defgroup interfaces_python Python interface
/// Interface to call Python functions.
/// @{

/// Polymorphic Python object
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
	object( const vec<T> &val );
	template<class T>
	object( const mat<T> &val );
	~object();
	void operator=( const object &o );
	void* data() const { return pdata; }
	operator int() const;
	operator double() const;
	operator complex() const;
	operator string() const;
	template<class T>
	operator vec<T>() const;
	template<class T>
	operator mat<T>() const;
};

/// Functor handling a Python function (for internal use)
class _function_handler {

	void* pFunc;																///< Pointer to the python function handler
	static int py_init_flag;													///< Flag to indicate if the interface has been initialized

	/// Initializes the interface.
	void init();

protected:

	/// Function call returning a polymorphic object
	object fcall( const deque<object> &O );

public:

	/// Constructor
	_function_handler( const string& module, const string& function );
};

/// Functor wrapper a Python function
template<class R>
class function : public _function_handler {
public:

	using _function_handler::_function_handler;

	/// Template function call
	template<class... T>
	R operator()( const T&... argin )
	{
	    // Parse parameters
	    deque<object> O;
	    O = parse_variadic_template<object>( argin... );

	    // Function call
	    R rv = fcall( O );
	    return rv;
	}
};

/// @}
}
}
#endif

#endif
