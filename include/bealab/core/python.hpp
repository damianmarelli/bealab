// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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
