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
		varlist.pop_back();
		varlist.pop_back();
		return varlist;
	}

protected:

	string call()
	{
		return "octave -q --eval";
	}

public:

	/// Constructor
	functor( const string& funname ) :
		mfilename(funname), matfilename(create_tmpfile()), mf( matfilename, true ) {}

	/// Destructor
	~functor() { std::remove( matfilename.data() ); }

	/// Constructor
	template<class... iargs>
	typename std::tuple_element<0,std::tuple<oargs...>>::type
	operator()( const iargs&... X )
	{
		// Parse input parameters
		string ivarlist = parse_input( X... );

		// Prepare output parameters
		int N = sizeof...(oargs);
		ostringstream osovarlist;
		for( int n = 0; n < N; n++ )
			osovarlist << "Y" << n << ", ";
		string ovarlist = osovarlist.str();
		if( N > 0 ) {
			ovarlist.pop_back();
			ovarlist.pop_back();
		}

		// Octave command
		ostringstream cmd;
		cmd << call() << " '"
			<< "load " << matfilename << "; "
			<< "[" << ovarlist << "] = " << mfilename << "( " << ivarlist << " ); "
			<< "save -v6 " << matfilename << " " << ovarlist << ";"
			<< "'";
//		cout << cmd.str() << endl;
		system(cmd.str());

		// Parse result
		typedef typename std::tuple_element<0,std::tuple<oargs...>>::type R;
		return mf.load<R>("Y0");
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
