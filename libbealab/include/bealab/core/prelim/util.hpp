/// @file bealab/core/prelim/util.hpp
/// A number of useful functions and classes.

#ifndef _BEALAB_PRELIM_UTIL_
#define	_BEALAB_PRELIM_UTIL_

#include <unistd.h>
#include <cxxabi.h>
#include <sys/time.h>

namespace bealab
{
/// @defgroup prelim_util Utilities
/// Some general purpose functions and classes.
/// @{

/// Macro to inherit constructors, until GCC implements it
#define INHERIT_CONSTRUCTORS( dervied, base )                \
    template<typename ...Args,                               \
             typename = typename std::enable_if              \
             <                                               \
                std::is_constructible<base, Args...>::value  \
             >::type>                                        \
    dervied( Args &&...args )                                \
        : base( std::forward<Args>(args)... ) {}             \

/// Macro to inherit assignment operators
#define INHERIT_ASSIGNMENTS( dervied, base )                 \
    template<typename ...Args,                               \
             typename = typename std::enable_if              \
             <                                               \
                std::is_assignable<base, Args...>::value     \
             >::type>                                        \
    dervied& operator=( const Args&... args )                \
    {                                                        \
		base::operator=( args... );                          \
		return *this;                                        \
	}                                                        \

/// Returns a string with the type of a variable.
template<class T>
string type_name( const T &x )
{
	return abi::__cxa_demangle(typeid(x).name(), 0, 0, NULL );
}

/// Prints an error message and aborts.
void error( const string &str );

/// Prints a warning message and continues.
void warning( const string &str );

/// System call managing the return value
int system( const string& call );

/// Create a pipe for reading and writing
int rw_pipe( const string& program, const vector<string>& arguments,
		int* write_fd, int* read_fd, int* error_fd=NULL );

/// @name Parse variadic templates

/// Internal function used by parse_variadic_template().
template<class S>
void _parse_variadic_template( deque<S>& V ) {}

/// Internal function used by parse_variadic_template().
template<class S, class R, class... T>
void _parse_variadic_template( deque<S>& V, const R& r, const T&... t )
{
	S s = r;
	V.push_back( s );
    _parse_variadic_template( V, t... );
}

/// Converts variadic template arguments into deque<S>.
template<class S, class... T>
deque<S> parse_variadic_template( const T&... t )
{
	deque<S> V;
	_parse_variadic_template( V, t... );
	return V;
}
/// @}

/// A timer to measure execution time.
class timer {
	double last_time;
public:
	timer()
	{
	    timeval time;
	    gettimeofday(&time, NULL);
	    last_time = time.tv_sec + time.tv_usec / 1e6;
	}
	double reset()
	{
	    timeval time;
	    gettimeofday(&time, NULL);
		double ctime = time.tv_sec + time.tv_usec / 1e6;
		double etime = ctime - last_time;
		last_time    = ctime;
		return etime;
	}
	double elapsed()
	{
	    timeval time;
	    gettimeofday(&time, NULL);
		double ctime = time.tv_sec + time.tv_usec / 1e6;
		double etime = ctime - last_time;
		return etime;
	}
};

/// @}
}
#endif
