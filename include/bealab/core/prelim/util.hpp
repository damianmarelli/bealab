// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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
