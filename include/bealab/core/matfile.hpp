// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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
	void* create_struct( const string& varname, const Vec<void*>& fields );		///< Create a MATIO structure
	static
	void* create_cell( const string& varname, const Mat<void*>& entries );		///< Create a MATIO cell array
	/// @}

	/// @name Parse a MATIO container
	static
	void* parse_struct( void* pvar, const string& field );						///< Parse a MATIO structure
	static
	Mat<void*> parse_cell( void* pvar );										///< Parse a MATIO sell array
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

/// Specialization of matvar for Vec<double>
template<>
struct matvar<Vec<double>> {
	static void* create( const string& varname, const Vec<double>& x );
	static Vec<double> parse( void* pvar );
};

/// Specialization of matvar for Vec<bool>
template<>
struct matvar<Vec<bool>> {
	static void* create( const string& varname, const Vec<bool>& x );
	static Vec<bool> parse( void* pvar );
};

/// Specialization of matvar for Vec<int>
template<>
struct matvar<Vec<int>> {
	static void* create( const string& varname, const Vec<int>& x );
	static Vec<int> parse( void* pvar );
};

/// Specialization of matvar for Vec<complex>
template<>
struct matvar<Vec<complex>> {
	static void* create( const string& varname, const Vec<complex>& x );
	static Vec<complex> parse( void* pvar );
};

/// Specialization of matvar for Mat<double>
template<>
struct matvar<Mat<double>> {
	static void* create( const string& varname, const Mat<double>& x );
	static Mat<double> parse( void* pvar );
};

/// Specialization of matvar for Mat<bool>
template<>
struct matvar<Mat<bool>> {
	static void* create( const string& varname, const Mat<bool>& x );
	static Mat<bool> parse( void* pvar );
};

/// Specialization of matvar for Mat<int>
template<>
struct matvar<Mat<int>> {
	static void* create( const string& varname, const Mat<int>& x );
	static Mat<int> parse( void* pvar );
};

/// Specialization of matvar for Mat<complex>
template<>
struct matvar<Mat<complex>> {
	static void* create( const string& varname, const Mat<complex>& x );
	static Mat<complex> parse( void* pvar );
};

/// Specialization of matvar for Vec<T>
template<class T>
struct matvar<Vec<T>> {
	static
	void* create( const string& varname, const Vec<T>& x )
	{
		int N = x.size();
		Mat<void*> entries(N,1);
		for( int n = 0; n < N; n++ )
			entries(n,0) = matvar<T>::create( "", x(n) );
		return matfile::create_cell( varname, entries );
	}

	static
	Vec<T> parse( void* pvar )
	{
		Mat<void*> entries = matfile::parse_cell( pvar );
		int N = entries.size1();
		Vec<T> x(N);
		for( int n = 0; n < N; n++ )
			x(n) = matvar<T>::parse( entries(n,0) );
		return x;
	}
};

/// Specialization of matvar for Mat<T>
template<class T>
struct matvar<Mat<T>> {
	static
	void* create( const string& varname, const Mat<T>& x )
	{
		int M = x.size1();
		int N = x.size2();
		Mat<void*> entries(M,N);
		for( int m = 0; m < M; m++ )
			for( int n = 0; n < N; n++ )
				entries(m,n) = matvar<T>::create( "", x(m,n) );
		return matfile::create_cell( varname, entries );
	}

	static
	Mat<T> parse( void* pvar )
	{
		Mat<void*> entries = matfile::parse_cell( pvar );
		int M = entries.size1();
		int N = entries.size2();
		Mat<T> x(M,N);
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
