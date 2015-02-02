/// @file bealab/core/matlab.hpp
/// Interface to call Matlab functions.

#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOMATLAB

#ifndef _BEALAB_MATLAB_
#define	_BEALAB_MATLAB_

#include <bealab/core/blas.hpp>

namespace bealab
{
namespace matlab
{

//==============================================================================
// INTERFASE
//==============================================================================
/// @defgroup interfaces_matlab Matlab interface
/// Interface to call Matlab functions.
/// @{

/// @name Internal functions
void _mlinit();						///< Initialize Matlab interface.
void _matlab_call( const char * );	///< Matlab call using a string.
/// @}

template<int> class functor;

/// Interpreter class.
class object {

	template<int>
	friend class functor;

	void *pdata;																///< Pointer to an mxArray
	bool own;																	///< Flag to indicate if the array is owned
	void create_array( const void*, bool, const ivec & );						///< Associate a new mxArray
	void fill_array( const void*, bool, const ivec &, ivec );					///< Fill an mxArray with data

public:
	void put( const char* ) const;												///< Put the mxArray into the Matlab sesion
	void get( const char* );													///< Get an mxArray from the Matlab sesion

	/// @name Constructors & destructor
	object() : pdata(NULL), own(false) {}										///< Empty object
	object( void *p ) : pdata(p), own(false) {}									///< Associate object with eixsting mxArray
	object( double*, ivec );
	object( complex*, ivec );
	object( const object& );
	template<class T>
	object( const T& x ) : own(false) { *this = x; }
	~object();
	/// @}

	/// @name Assignment
	template<class T>
	object& operator=( T );
	template<class T>
	object& operator=( const vec<T>& );
	template<class T>
	object& operator=( const mat<T>& );
	object& operator=( const object& );
	object& operator=( const char* );
	object& operator=( const string & );
	/// @}

	/// @name Cast operators
//	template<class T,
//		class=typename enable_if<
//			is_fundamental<T>::value || is_same<T,complex>::value>
//		::type>
	template<class T>
	operator T() const;
	template<class T>
	operator vec<T>() const;
	template<class T>
	operator mat<T>() const;
	/// @}

	/// @name Pointer manipulation
	void* data() const;
	void data(void*);
	/// @}
};

///// Class for choosing the return type of the Matlab functor wrapper
//template<int N>
//struct choose_t {
//	typedef deque<object> type;
//};
//
///// Specialization for no return values
//template<>
//struct choose_t<0> {
//	typedef void type;
//};
//
///// Specialization for one return value
//template<>
//struct choose_t<1> {
//	typedef object type;
//};

/// Functor wrapper for a Matlab function
template<int noutdefault=-1>
class functor {

	/// Name of the Matlab function
	string funname;

//	/// Function call with a specified number of output arguments
//	template<class... T>
//	deque<object> fcall( const T&... argin )
//	{
//		// Parsing
//		deque<object> varg;
//		varg = parse_variadic_template<object>( argin... );
//
//		// Start engine
//	    _mlinit();
//
//		// Form command line ------------------------------------------------
//		string fcall;
//		if( nargout > 0) {
//			fcall += "[";
//
//			// Output arguments
//			for( int j = 0; j < nargout; j++ ) {
//
//				// variable name
//				std::ostringstream varname;
//				varname << "y" << j;
//
//				// form function call
//				fcall += varname.str() + ", ";
//			}
//			fcall.erase( fcall.length()-2 );
//			fcall += "] = ";
//		}
//
//		//function call
//		 fcall += funname + "( ";
//
//		// Input arguments
//		int I = varg.size();
//		for (int i = 0; i < I; i++) {
//
//			// variable name
//			std::ostringstream varname;
//			varname << "x" << i;
//
//			// put variable in Matlab
//			varg[i].put( varname.str().data() );
//
//			// form function call
//			fcall += varname.str() + ", ";
//		}
//		if( I > 0 )
//		fcall.erase( fcall.length()-2 );
//
//		fcall += " );";
//
//		// Matlab function call --------------------------------------------
//		_matlab_call( fcall.data() );
//
//		// Retrieve command output -----------------------------------------
//		deque<object> argout;
//		for( int j = 0; j < nargout; j++ ) {
//
//			// variable name
//			std::ostringstream varname;
//			varname << "y" << j;
//
//			// get variable from Matlab
//			object y;
//			y.get( varname.str().data() );
//			argout.push_back( y );
//		}
//
//		return argout;
//	}

	int get_nargout()
	{
		if( noutdefault != -1 )
			return noutdefault;
		else {
			string fcall = "n = nargout('"+funname+"');";
			_matlab_call( fcall.data() );
			object n;
			n.get( "n" );
			int in = n;
			if( in < 0 )
				error("matlab::functor - the function '"+funname+"' has a variable number of arguments, so a template parameter needs to be specified");
			return in;
		}
	}

public:

	/// Constructor
	functor( const string &fn ) : funname(fn) {};

	/// Call operator
	template<class... T>
	deque<object> operator()( const T&... argin )
	{
		// Parsing
		deque<object> varg;
		varg = parse_variadic_template<object>( argin... );

		// Start engine
		_mlinit();

		// Form command line ------------------------------------------------
		// Output arguments
		int nargout = get_nargout();
		string fcall;
		if( nargout > 0) {
			fcall += "[";

			// Output arguments
			for( int j = 0; j < nargout; j++ ) {

				// variable name
				std::ostringstream varname;
				varname << "y" << j;

				// form function call
				fcall += varname.str() + ", ";
			}
			fcall.erase( fcall.length()-2 );
			fcall += "] = ";
		}

		//function call
		 fcall += funname + "( ";

		// Input arguments
		int I = varg.size();
		for (int i = 0; i < I; i++) {

			// variable name
			std::ostringstream varname;
			varname << "x" << i;

			// put variable in Matlab
			varg[i].put( varname.str().data() );

			// form function call
			fcall += varname.str() + ", ";
		}
		if( I > 0 )
		fcall.erase( fcall.length()-2 );

		fcall += " );";

		// Matlab function call --------------------------------------------
		_matlab_call( fcall.data() );

		// Retrieve command output -----------------------------------------
		deque<object> argout;
		for( int j = 0; j < nargout; j++ ) {

			// variable name
			std::ostringstream varname;
			varname << "y" << j;

			// get variable from Matlab
			object y;
			y.get( varname.str().data() );
			argout.push_back( y );
		}

		return argout;
	}
};

///// Specialization for no return values
//template<>
//template<class... T>
//typename choose_t<0>::type functor<0>::operator()( const T&... argin )
//{
//	fcall( argin... );
//}
//
///// Specialization for one return value
//template<>
//template<class... T>
//typename choose_t<1>::type functor<1>::operator()( const T&... argin )
//{
//	return fcall( argin... )[0];
//}

/// @name Interfase functions

/// Turns interface tracing on/off.
void trace( bool );

/// @}

/// @name Plotting interface
static functor<> figure( "figure" );
static functor<> close( "close" );
static functor<> plot( "plot" );
static functor<> semilogx( "semilogx" );
static functor<> semilogy( "semilogy" );
static functor<> loglog( "loglog" );
static functor<> mesh( "mesh" );
static functor<> hold( "hold" );
static functor<> title( "title" );
static functor<> axes( "axes" );
static functor<> axis( "axis" );
static functor<> xlabel( "xlabel" );
static functor<> ylabel( "ylabel" );
static functor<> xlim( "xlim" );
static functor<> ylim( "ylim" );
static functor<> subplot( "subplot" );
static functor<> legend( "legend" );
static functor<> imagesc( "imagesc" );
static functor<> colormap( "colormap" );
/// @}

/// @}
}
}
#endif

#endif

/*
 * TI VOGLIO TANTO BENE :* (smak!)
 */
