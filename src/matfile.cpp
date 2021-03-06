/*******************************************************************************
 * This software is licensed under the BSD 3-Clause License with the possibility to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
 * You may not use this work except in compliance with the License.
 * You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
 * Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the file License for the specific language governing permissions and limitations under the License. 
 * If you wish to obtain a commercial license, please contact the authors via e-mail.
 *
 * Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)
 *******************************************************************************/
#include <bealab/core/matfile.hpp>
#include <matio.h>
#include <dirent.h>
#include <fnmatch.h>

namespace bealab
{

//------------------------------------------------------------------------------
// matfile class definition
//------------------------------------------------------------------------------
void matfile::_open()
{
    pfile = Mat_Open( filename.data(), MAT_ACC_RDWR );
    if( pfile == NULL )
		error( "matfile::_open() - Error opening file " + filename );
}

void matfile::_close()
{
	Mat_Close( (mat_t*)pfile );
}

void* matfile::_read( const string& varname )
{
	void* rv = Mat_VarRead( (mat_t*)pfile, const_cast<char*>(varname.data()) );
	if( rv == NULL ) {
		if( trace == true )
			cerr << "Variable '" << varname << "' not found in file '"
			     << filename << "'." << endl;
		throw variable_not_found(varname);
	}
	return rv;
}

void matfile::_write( void* pmatvar )
{
	if( Mat_VarWrite ( (mat_t*)pfile, (matvar_t*)pmatvar, MAT_COMPRESSION_NONE ) != 0 )
		error("matfile::_write() - Error writing variable");
}

// XXX This function needs to free the members of a struct/cell
void matfile::_free( void* pmatvar )
{
//	matvar_t* pvar   = reinterpret_cast<matvar_t*>(pmatvar);
//	matvar_t** pdata = reinterpret_cast<matvar_t**>( pvar->data );
//	switch( pvar->class_type ) {
//	case MAT_C_STRUCT:
//		for( int i = 0;; i++ )
//			if( pdata[i] != NULL )
//				_free( pdata[i] );
//			else
//				break;
//		break;
//	case MAT_C_CELL:
//		int I = pvar->dims[0] * pvar->dims[1];
//		for( int i = 0; i < I; i++ )
//			_free( pdata[i] );
//		break;
//	}
//	Mat_VarFree( pvar );
	Mat_VarFree( (matvar_t*)pmatvar );
}

void* matfile::create_struct( const string& varname, const vec<void*>& fields )
{
	size_t dims[2] = {1,1};
	int N = fields.size();
	matvar_t *matvar[N+1];
	for( int n = 0; n < N; n++ )
		matvar[n] = (matvar_t*)fields(n);
	matvar[N] = NULL;
	return Mat_VarCreate( varname.data(), MAT_C_STRUCT, MAT_T_STRUCT, 2, dims, matvar, 0 );
}

void* matfile::create_cell( const string& varname, const mat<void*>& entries )
{
	size_t I = entries.size1();
	size_t J = entries.size2();
	size_t dims[2] = {I,J};
	matvar_t *matvar[J][I];
	for( size_t i = 0; i < I; i++ )
		for( size_t j = 0; j < J; j++ )
			matvar[j][i] = (matvar_t*)entries(i,j);
	return Mat_VarCreate( varname.data(), MAT_C_CELL, MAT_T_CELL, 2, dims, matvar, 0 );
}

void* matfile::parse_struct( void* pvar, const string& field )
{
	return Mat_VarGetStructField( (matvar_t*)pvar, const_cast<char*>(field.data()), MAT_BY_NAME, 0 );
}

mat<void*> matfile::parse_cell( void* pvar )
{
	matvar_t* pmatvar = (matvar_t*)pvar;
	int I = pmatvar->dims[0];
	int J = pmatvar->dims[1];
	mat<void*> rv(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			rv(i,j) = reinterpret_cast<matvar_t**>(pmatvar->data)[i+j*I];
	return rv;
}

matfile::matfile( const string& filename_, bool create ) : filename(filename_)
{
	if( create || !bealab::exists(filename) )
        pfile = Mat_Create( filename.data(), NULL );
	else
		pfile = Mat_Open( filename.data(), MAT_ACC_RDWR );
	if( pfile == NULL )
        	error( "matfile::matfile() - Error creating/opening file " + filename );
	Mat_Close( (mat_t*)pfile );
}

bool matfile::exists( const string& varname )
{
	_open();
	void* rv = Mat_VarRead( (mat_t*)pfile, const_cast<char*>(varname.data()) );
	_close();
	if( rv == NULL )
		return false;
	else
		return true;
}

void matfile::erase( const string& varname )
{
	_open();
	Mat_VarDelete ( (mat_t*)pfile, const_cast<char*>(varname.data()) );
	_close();
}

//------------------------------------------------------------------------------
// matvar<string> struct definition
//------------------------------------------------------------------------------
void* matvar<string>::create( const string& varname, const string& x )
{
	// Create variable
	size_t I       = x.size();
	size_t dims[2] = {1,I};
	char data[I];
	for( size_t i = 0; i < I; i++ )
		data[i] = x.data()[i];
	matvar_t* var = Mat_VarCreate( varname.data(), MAT_C_CHAR, MAT_T_UINT8,
								   2, dims, data, 0 );
	// Return the pointer
	return var;
}

string matvar<string>::parse( void* pvar )
{
	matvar_t* pmatvar = (matvar_t*)pvar;
	auto pdata        = reinterpret_cast<char*>(pmatvar->data);
	return pdata;
}

//------------------------------------------------------------------------------
// matvar<double> struct definition
//------------------------------------------------------------------------------
void* matvar<double>::create( const string& varname, const double& x )
{
	// Create variable
	size_t dims[2] = {1,1};
	double data[1] = {x};
	matvar_t* var = Mat_VarCreate( varname.data(), MAT_C_DOUBLE, MAT_T_DOUBLE,
								   2, dims, data, 0 );
	// Return the pointer
	return var;
}

double matvar<double>::parse( void* pvar )
{
	matvar_t* pmatvar = (matvar_t*)pvar;
	auto pdata    = reinterpret_cast<double*>(pmatvar->data);
	return pdata[0];
}

//------------------------------------------------------------------------------
// matvar<bool> struct definition
//------------------------------------------------------------------------------
void* matvar<bool>::create( const string& varname, const bool& x )
{
	return matvar<double>::create( varname, x );
}

bool matvar<bool>::parse( void* pvar )
{
	return matvar<double>::parse( pvar );
}

//------------------------------------------------------------------------------
// matvar<int> struct definition
//------------------------------------------------------------------------------
void* matvar<int>::create( const string& varname, const int& x )
{
	return matvar<double>::create( varname, x );
}

int matvar<int>::parse( void* pvar )
{
	return matvar<double>::parse( pvar );
}

//------------------------------------------------------------------------------
// matvar<complex> struct definition
//------------------------------------------------------------------------------
void* matvar<complex>::create( const string& varname, const complex& x )
{
	// Create variable
	size_t dims[2] = {1,1};
	double r = real(x);
	double i = imag(x);
	mat_complex_split_t z = {&r,&i};
	matvar_t* var = Mat_VarCreate( varname.data(), MAT_C_DOUBLE, MAT_T_DOUBLE,
								   2, dims, &z, MAT_F_COMPLEX );
	// Return the pointer
	return var;
}

complex matvar<complex>::parse( void* pvar )
{
	matvar_t* pmatvar = (matvar_t*)pvar;
	auto pdata    = reinterpret_cast<mat_complex_split_t*>(pmatvar->data);
	auto preal    = reinterpret_cast<double*>(pdata->Re);
	auto pimag    = reinterpret_cast<double*>(pdata->Im);
	return complex{ preal[0], pimag[0] };
}

//------------------------------------------------------------------------------
// matvar<vec<double>> struct definition
//------------------------------------------------------------------------------
void* matvar<vec<double>>::create( const string& varname, const vec<double>& x )
{
	// Create variable
	size_t I = x.size();
	size_t dims[2] = {I,1};
	double data[I];
	for( size_t i = 0; i < I; i++ )
		data[i] = x(i);
	matvar_t* var = Mat_VarCreate( varname.data(), MAT_C_DOUBLE, MAT_T_DOUBLE,
								   2, dims, data, 0 );
	// Return the pointer
	return var;
}

vec<double> matvar<vec<double>>::parse( void* pvar )
{
	matvar_t* pmatvar = (matvar_t*)pvar;
	auto pdata    = reinterpret_cast<double*>(pmatvar->data);
	int I = pmatvar->dims[0];
	vec<double> x(I);
	for( int i = 0; i < I; i++ )
		x(i) = pdata[i];
	return x;
}

//------------------------------------------------------------------------------
// matvar<vec<bool>> struct definition
//------------------------------------------------------------------------------
void* matvar<vec<bool>>::create( const string& varname, const vec<bool>& x )
{
	return matvar<vec<double>>::create( varname, x );
}

vec<bool> matvar<vec<bool>>::parse( void* pvar )
{
	return matvar<vec<double>>::parse( pvar );
}

//------------------------------------------------------------------------------
// matvar<vec<int>> struct definition
//------------------------------------------------------------------------------
void* matvar<vec<int>>::create( const string& varname, const vec<int>& x )
{
	return matvar<vec<double>>::create( varname, x );
}

vec<int> matvar<vec<int>>::parse( void* pvar )
{
	return matvar<vec<double>>::parse( pvar );
}

//------------------------------------------------------------------------------
// matvar<vec<complex>> struct definition
//------------------------------------------------------------------------------
void* matvar<vec<complex>>::create( const string& varname, const vec<complex>& x )
{
	// Create variable
	size_t I = x.size();
	size_t dims[2] = {I,1};
	double datar[I];
	double datai[I];
	for( size_t i = 0; i < I; i++ ) {
		datar[i] = real(x(i));
		datai[i] = imag(x(i));
	}
	struct mat_complex_split_t z = {datar,datai};
	matvar_t* var = Mat_VarCreate( varname.data(), MAT_C_DOUBLE, MAT_T_DOUBLE,
								   2, dims, &z, (I == 0 ? 0 : MAT_F_COMPLEX) );
	// Return the pointer
	return var;
}

vec<complex> matvar<vec<complex>>::parse( void* pvar )
{
	matvar_t* pmatvar = (matvar_t*)pvar;
	auto pdata    = reinterpret_cast<mat_complex_split_t*>(pmatvar->data);
	int I = pmatvar->dims[0];
	vec<complex> x(I);
	for( int i = 0; i < I; i++ ) {
		auto preal = reinterpret_cast<double*>(pdata->Re);
		auto pimag = reinterpret_cast<double*>(pdata->Im);
		x(i)       = preal[i] + bealab::i * pimag[i];
	}
	return x;
}

//------------------------------------------------------------------------------
// matvar<mat<double>> struct definition
//------------------------------------------------------------------------------
void* matvar<mat<double>>::create( const string& varname, const mat<double>& x )
{
	// Create variable
	size_t I = x.size1();
	size_t J = x.size2();
	size_t dims[2] = {I,J};
	double data[J][I];
	for( size_t i = 0; i < I; i++ )
		for( size_t j = 0; j < J; j++ )
			data[j][i] = x(i,j);
	matvar_t* var = Mat_VarCreate( varname.data(), MAT_C_DOUBLE, MAT_T_DOUBLE,
								   2, dims, data, 0 );
	// Return the pointer
	return var;
}

mat<double> matvar<mat<double>>::parse( void* pvar )
{
	matvar_t* pmatvar = (matvar_t*)pvar;
	auto pdata    = reinterpret_cast<double*>(pmatvar->data);
	int I = pmatvar->dims[0];
	int J = pmatvar->dims[1];
	mat<double> x(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			x(i,j) = pdata[i+j*I];
	return x;
}

//------------------------------------------------------------------------------
// matvar<mat<bool>> struct definition
//------------------------------------------------------------------------------
void* matvar<mat<bool>>::create( const string& varname, const mat<bool>& x )
{
	return matvar<mat<double>>::create( varname, x );
}

mat<bool> matvar<mat<bool>>::parse( void* pvar )
{
	return matvar<mat<double>>::parse( pvar );
}

//------------------------------------------------------------------------------
// matvar<mat<int>> struct definition
//------------------------------------------------------------------------------
void* matvar<mat<int>>::create( const string& varname, const mat<int>& x )
{
	return matvar<mat<double>>::create( varname, x );
}

mat<int> matvar<mat<int>>::parse( void* pvar )
{
	return matvar<mat<double>>::parse( pvar );
}

//------------------------------------------------------------------------------
// matvar<mat<complex>> struct definition
//------------------------------------------------------------------------------
void* matvar<mat<complex>>::create( const string& varname, const mat<complex>& x )
{
	// Create variable
	size_t I = x.size1();
	size_t J = x.size2();
	size_t dims[2] = {I,J};
	double datar[J][I];
	double datai[J][I];
	for( size_t i = 0; i < I; i++ )
		for( size_t j = 0; j < J; j++ ) {
			datar[j][i] = real(x(i,j));
			datai[j][i] = imag(x(i,j));
		}
	mat_complex_split_t z = {datar,datai};
	matvar_t* var = Mat_VarCreate( varname.data(), MAT_C_DOUBLE, MAT_T_DOUBLE,
								   2, dims, &z, ((I*J) == 0 ? 0 : MAT_F_COMPLEX) );
	// Return the pointer
	return var;
}

mat<complex> matvar<mat<complex>>::parse( void* pvar )
{
	matvar_t* pmatvar = (matvar_t*)pvar;
	auto pdata    = reinterpret_cast<mat_complex_split_t*>(pmatvar->data);
	int I = pmatvar->dims[0];
	int J = pmatvar->dims[1];
	mat<complex> x(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ ) {
			auto preal = reinterpret_cast<double*>(pdata->Re);
			auto pimag = reinterpret_cast<double*>(pdata->Im);
			x(i,j)     = preal[i+j*I] + bealab::i * pimag[i+j*I];
		}
	return x;
}

bool exists( const string& filename )
{
	FILE* fp = fopen( filename.data(), "r" );
	if( fp == NULL )
		return false;
	fclose( fp );
	return true;
}

deque<string> listdir( const string& dirname, const string& pattern )
{
	deque<string> fnames;
	DIR* dir = opendir ( dirname.data() );
	struct dirent *entry;
	while( (entry = readdir(dir)) )
		if( !fnmatch( pattern.data(), entry->d_name, FNM_CASEFOLD ) )
			fnames.push_back( dirname + entry->d_name );
	closedir (dir);
	return fnames;
}

}
