/// @file bealab/core/matfile.hpp
/// Interface for loading and saving bealab objects in MAT format

#ifndef _BEALAB_MATFILE_
#define	_BEALAB_MATFILE_

#include <exception>
#include <bealab/core/blas.hpp>

namespace bealab
{
/// @defgroup matfile Load and save
/// Interface for loading and saving bealab objects in MAT format
/// @{

/// Interpreter structure for a matvar
template<class T> struct matvar;

/// Class for load/save to from a matfile
class matfile {

	string filename;															///< File name.
	void* pfile;																/// Pointer to the matfile handler

	/// @name Low level matfile access
	void _open();																///< Open the matfile to put or read a variable
	void _close();																///< Close the matfile after read/write
	void* _read( const string& varname );										///< Read a variable from the matfile
	void _write( void* pmatvar );												///< Write a variable to the matfile
	/// @}

	/// @name Free a MATIO variable
	void _free( void* pmatvar );												///< Destroy a MATIO variable
	/// @}

public:

	bool trace = true;															///< Trace when a variable is not found

	/// Exception thrown when a variable is not found
	class variable_not_found : public std::exception {
	public:
		string varname;
		variable_not_found( const string& s ) : varname(s) {}
		~variable_not_found() throw() {}
	};

	/// @name Create a MATIO container
	static
	void* create_struct( const string& varname, const vec<void*>& fields );		///< Create a MATIO structure
	static
	void* create_cell( const string& varname, const mat<void*>& entries );		///< Create a MATIO cell array
	/// @}

	/// @name Parse a MATIO container
	static
	void* parse_struct( void* pvar, const string& field );						///< Parse a MATIO structure
	static
	mat<void*> parse_cell( void* pvar );										///< Parse a MATIO sell array
	/// @}

	/// Constructor
	matfile( const string& filename_, bool create=false );

	/// Save a variable
	template<class T>
	matfile& save( const string& varname, const T& x )
	{
		// If the variable exists, delete it
		if( exists(varname) )
			erase( varname );

		// Create and save the variable
		void* pvar = matvar<T>::create( varname, x );
		_open();
		_write( pvar );
		_close();
		_free( pvar );

		return *this;
	}

	/// Load a variable
	template<class T>
	matfile& load( const string& varname, T& x )
	{
		_open();
		void* pvar;
		pvar = _read( varname );
		_close();
		x = matvar<T>::parse( pvar );
		_free( pvar );
		return *this;
	}

	/// Load a variable
	template<class T>
	T load( const string& varname )
	{
		T x;
		load( varname, x );
		return x;
	}

	/// Check if a variable exists in the matfile
	bool exists( const string& varname );

	/// Delete a variable from the matfile
	void erase( const string& varname );
};

/// @name matvar specializations

/// Specialization of matvar for strings
template<>
struct matvar<string> {
	static void* create( const string& varname, const string& x );
	static string parse( void* pvar );
};

/// Specialization of matvar for bool
template<>
struct matvar<bool> {
	static void* create( const string& varname, const bool& x );
	static bool parse( void* pvar );
};

/// Specialization of matvar for int
template<>
struct matvar<int> {
	static void* create( const string& varname, const int& x );
	static int parse( void* pvar );
};

/// Specialization of matvar for double
template<>
struct matvar<double> {
	static void* create( const string& varname, const double& x );
	static double parse( void* pvar );
};

/// Specialization of matvar for complex
template<>
struct matvar<complex> {
	static void* create( const string& varname, const complex& x );
	static complex parse( void* pvar );
};

/// Specialization of matvar for vec<double>
template<>
struct matvar<vec<double>> {
	static void* create( const string& varname, const vec<double>& x );
	static vec<double> parse( void* pvar );
};

/// Specialization of matvar for vec<bool>
template<>
struct matvar<vec<bool>> {
	static void* create( const string& varname, const vec<bool>& x );
	static vec<bool> parse( void* pvar );
};

/// Specialization of matvar for vec<int>
template<>
struct matvar<vec<int>> {
	static void* create( const string& varname, const vec<int>& x );
	static vec<int> parse( void* pvar );
};

/// Specialization of matvar for vec<complex>
template<>
struct matvar<vec<complex>> {
	static void* create( const string& varname, const vec<complex>& x );
	static vec<complex> parse( void* pvar );
};

/// Specialization of matvar for mat<double>
template<>
struct matvar<mat<double>> {
	static void* create( const string& varname, const mat<double>& x );
	static mat<double> parse( void* pvar );
};

/// Specialization of matvar for mat<bool>
template<>
struct matvar<mat<bool>> {
	static void* create( const string& varname, const mat<bool>& x );
	static mat<bool> parse( void* pvar );
};

/// Specialization of matvar for mat<int>
template<>
struct matvar<mat<int>> {
	static void* create( const string& varname, const mat<int>& x );
	static mat<int> parse( void* pvar );
};

/// Specialization of matvar for mat<complex>
template<>
struct matvar<mat<complex>> {
	static void* create( const string& varname, const mat<complex>& x );
	static mat<complex> parse( void* pvar );
};

/// Specialization of matvar for vec<T>
template<class T>
struct matvar<vec<T>> {
	static
	void* create( const string& varname, const vec<T>& x )
	{
		int N = x.size();
		mat<void*> entries(N,1);
		for( int n = 0; n < N; n++ )
			entries(n,0) = matvar<T>::create( "", x(n) );
		return matfile::create_cell( varname, entries );
	}

	static
	vec<T> parse( void* pvar )
	{
		mat<void*> entries = matfile::parse_cell( pvar );
		int N = entries.size1();
		vec<T> x(N);
		for( int n = 0; n < N; n++ )
			x(n) = matvar<T>::parse( entries(n,0) );
		return x;
	}
};

/// Specialization of matvar for mat<T>
template<class T>
struct matvar<mat<T>> {
	static
	void* create( const string& varname, const mat<T>& x )
	{
		int M = x.size1();
		int N = x.size2();
		mat<void*> entries(M,N);
		for( int m = 0; m < M; m++ )
			for( int n = 0; n < N; n++ )
				entries(m,n) = matvar<T>::create( "", x(m,n) );
		return matfile::create_cell( varname, entries );
	}

	static
	mat<T> parse( void* pvar )
	{
		mat<void*> entries = matfile::parse_cell( pvar );
		int M = entries.size1();
		int N = entries.size2();
		mat<T> x(M,N);
		for( int m = 0; m < M; m++ )
			for( int n = 0; n < N; n++ )
				x(m,n) = matvar<T>::parse( entries(m,n) );
		return x;
	}
};
/// @}

/// Check if a file exists
bool exists( const string& filename );

/// List directory entries matching a pattern
deque<string> listdir( const string& dirname, const string& pattern="*" );

/// @}
}
#endif
