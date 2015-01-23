// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/core/octave.hpp
/// Interface to call Octave functions.

#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOOCTAVE

#ifndef _BEALAB_OCTAVE_
#define	_BEALAB_OCTAVE_

#include <bealab/core/blas.hpp>
#include <bealab/core/matfile.hpp>

namespace bealab
{
/// Octave interface
namespace octave
{
/// @defgroup interfaces_octave Octave interface
/// Interface to call Octave functions.
/// @{

/// Functor holding an Octave function
template<class... oargs>
class functor {

	string mfilename;															///< Name of the *.m file to containing the Octave function
	string tmpfilename;															///< Name of the associated temporary file
	string matfilename;															///< Name of the *.mat file to exchange data
	matfile mf;																	///< Handler of the *.mat file

	/// Create a temporary filename
	string create_tmpfile()
	{
		char fn[] = "/tmp/matfileXXXXXX";
		int fd = mkstemp( fn );
		if( fd == -1 )
			error("Cannot create temporary file");
		if( close(fd) == -1 )
			error("Cannot close temporary file");
		return fn;
	}

	void parse_input_( int n, string& varlist ) {}

	template<class T, class... S>
	void parse_input_( int n, string& varlist, const T& x, const S&... X )
	{
		ostringstream varname;
		varname << "X" << n;
		mf.save( varname.str(), x );
		parse_input_( n+1, varlist, X... );
		varlist = varname.str() + ", " + varlist;
	}

	template<class... T>
	string parse_input( const T&... X )
	{
		string varlist;
		parse_input_( 0, varlist, X... );
		if( sizeof...(X) > 0 ) {
			varlist.pop_back();
			varlist.pop_back();
		}
		return varlist;
	}

	template<class R = void>
	typename enable_if<std::is_void<R>::value,tuple<>>::type
	parse_output_( int n )
	{
		return tuple<>();
	}

	template<class T, class... S>
	tuple<T,S...> parse_output_( int n )
	{
		ostringstream varname;
		varname << "Y" << n;
		T x           = mf.load<T>(varname.str());
		tuple<S...> X = parse_output_<S...>( n+1 );
		return tuple_cat( tuple<T>(x), X );
	}

	tuple<oargs...> parse_output()
	{
		return parse_output_<oargs...>( 0 );
	}

protected:

	virtual
	string call_prefix()
	{
		return "octave -q --eval";
	}

	virtual
	string call_suffix()
	{
		return string();
	}

public:

	/// Constructor
	functor( const string& funname ) :
		mfilename(funname),
		tmpfilename(create_tmpfile()),
		matfilename(tmpfilename + ".mat"),
		mf( matfilename, true )
	{}

	/// Destructor
	~functor()
	{
		std::remove( tmpfilename.data() );
		std::remove( matfilename.data() );
	}

	/// Constructor
	template<class... iargs>
	tuple<oargs...> operator()( const iargs&... X )
	{
		// Parse input parameters
		string ivarlist = parse_input( X... );

		// Prepare output parameters
		int N = sizeof...(oargs);
		ostringstream osovarlist, ossvarlist;
		for( int n = 0; n < N; n++ ) {
			osovarlist << "Y" << n << ", ";
			ossvarlist << "Y" << n << " ";
		}
		string ovarlist = osovarlist.str();
		string svarlist = ossvarlist.str();
		if( N > 0 ) {
			ovarlist.pop_back();
			ovarlist.pop_back();
		}

		// Octave command
		ostringstream cmd;
		cmd << call_prefix() << " '"
			<< "load " << matfilename << "; " ;
		if( N > 0 )
			cmd<< "[" << ovarlist << "] = ";
		cmd << mfilename << "( " << ivarlist << " ); ";
		if( N > 0 )
			cmd << "save -v6 " << matfilename << " " << svarlist << ";";
		cmd << "exit;' " << call_suffix();
//		cout << cmd.str() << endl;
		system(cmd.str());

		// Parse result
		return parse_output();
	}
};

/// @}
}
}
#endif
#endif

/*
 * TI VOGLIO TANTO BENE :* (smak!)
 */
