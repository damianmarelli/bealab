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
#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOMATLAB

#include <bealab/core/matlab.hpp>
#include <engine.h>

using namespace bealab;

//------------------------------------------------------------------------------
// Private
//------------------------------------------------------------------------------
static Engine *_mleng = NULL;

static bool _trace = false;

static
class init{
public:
	~init() { if( _mleng != NULL ) engClose( _mleng ); }
} _init;

//------------------------------------------------------------------------------
// Namespace for matlab
//------------------------------------------------------------------------------
namespace bealab
{
namespace matlab
{

//---------------------------------------------------------------------------
// Internal functions
//---------------------------------------------------------------------------
void _mlinit()
{
	if( _mleng == NULL ){
		if ( !(_mleng = engOpen("matlab -nosplash")) ) {
			std::cerr << "\nCan't start MATLAB engine. Make sure that 'csh' is installed." << std::endl;
			abort();
		}
	}
}

void _matlab_call( const char *command )
{
	int N = 10e3;
	char buffer[N+1];
	buffer[N] = '\0';
	engOutputBuffer( _mleng, buffer, N );
	if( _trace )
		cout << command << endl;
	engEvalString( _mleng, command );
	if( _trace )
		cout << buffer << endl;
}

//---------------------------------------------------------------------------
// object class
//---------------------------------------------------------------------------
void object::fill_array( const void *psource, bool cmp, const ivec &dimsleft, ivec idx )
{
	int I    = dimsleft( 0 );
	int ndim = idx.size() - dimsleft.size();

	// If still more than one dimension, reduce one
	if( dimsleft.size() > 1 ) {
		ivec rdims = dimsleft( range( 1, dimsleft.size() ) );
		int skip   = prod( rdims );
		for( int i = 0; i < I; i++ ) {
			idx(ndim) = i;
			if( cmp == false )
				fill_array( static_cast<const double*>(psource)+i*skip, cmp, rdims, idx );
			else
				fill_array( static_cast<const complex*>(psource)+i*skip, cmp, rdims, idx );
		}
	}
	// Otherwise convert the data
	else {
		double *preal = mxGetPr( (mxArray*)pdata );
		double *pimag = NULL;
		if( cmp == true )
			pimag = mxGetPi( (mxArray*)pdata );
		vec<mwIndex> idx_ = idx;
		for( int i = 0; i < I; i++ ) {
			idx_(ndim) = i;
			int offset = mxCalcSingleSubscript( (mxArray*)pdata, idx_.size(), idx_.data().begin() );
			if( cmp == false )
				preal[offset] = static_cast<const double*>(psource)[i];
			else {
				preal[offset] = real(static_cast<const complex*>(psource)[i]);
				pimag[offset] = imag(static_cast<const complex*>(psource)[i]);
			}
		}
	}
}

void object::create_array( const void *psource, bool cmp, const ivec &dims )
{
	// Destroy a possible previous array
	if( own )
		mxDestroyArray( (mxArray*)pdata );

	// Create a new mxArray
	vec<mwSize> dims_ = dims;
	mxComplexity cmpl;
	if( cmp == false )
		cmpl = mxREAL;
	else
		cmpl = mxCOMPLEX;
	pdata = mxCreateNumericArray( dims_.size(), dims_.data().begin(),
							mxDOUBLE_CLASS, cmpl );

	// Fill it with the data pointed by pdata
	ivec idx(dims.size());
	fill_array( psource, cmp, dims, idx );

	// Mark the mxArray as owned
	own = true;
}

// Put in matlab
void object::put( const char *name ) const
{
	_mlinit();
	engPutVariable( _mleng, name, (mxArray*)pdata );
}

// Get from matlab
void object::get( const char *name )
{
	// Destroy a possible previous array
	if( own )
		mxDestroyArray( (mxArray*)pdata );

	// Get a new mxArray from Matlab
	_mlinit();
	pdata = engGetVariable( _mleng, name );
	if( pdata == NULL )
		error("Error using matGetVariable()");

	// Mark the mxArray as not owned
	own = true;
}

// Constructors & destructor ---------------------------------------------------
object::object( const object& x ) : own(false)
{
	*this = x;
}

object::object( double *psource, ivec dims )
{
	create_array( psource, false, dims );
}

object::object( complex *psource, ivec dims )
{
	create_array( psource, true, dims );
}

object::~object()
{
	// Destroy a possible previous array
	if( own )
		mxDestroyArray( (mxArray*)pdata );
}

// Assignment ------------------------------------------------------------------
template<class T>
object& object::operator=( T v_ )
{
	double v = v_;
	ivec dims = { 1 };
	create_array( &v, false, dims );
	return *this;
}

template<class T>
object& object::operator=( const vec<T> &v_ )
{
	rvec v = v_;
	ivec dims = { (int)v.size() };
	create_array( v.data().begin(), false, dims );
	return *this;
}

template<class T>
object& object::operator=( const mat<T> &m_ )
{
	rmat m = m_;
	ivec dims = { (int)m.size1(), (int)m.size2() };
	create_array( m.data().begin(), false, dims );
	return *this;
}

template<>
object& object::operator=( complex v_ )
{
	complex v = v_;
	ivec dims = { 1 };
	create_array( &v, true, dims );
	return *this;
}

template<>
object& object::operator=( const cvec &v_ )
{
	cvec v = v_;
	ivec dims = { (int)v.size() };
	create_array( v.data().begin(), true, dims );
	return *this;
}

template<>
object& object::operator=( const cmat &m_ )
{
	cmat m = m_;
	ivec dims = { (int)m.size1(), (int)m.size2() };
	create_array( m.data().begin(), true, dims );
	return *this;
}

object& object::operator=( const object &a )
{
	// Destroy a possible previous array
	if( own )
		mxDestroyArray( (mxArray*)pdata );

	// Copy the whole mxArray and mark it as owned
	pdata= mxDuplicateArray( (mxArray*)a.pdata );
	own = true;

//	a.put( "tmp" );
//	get( "tmp" );
	return *this;
}

object& object::operator=( const char *str )
{
	// Destroy a possible previous array
	if( own )
		mxDestroyArray( (mxArray*)pdata );

	// Make a new object an dmark it as owned
	pdata = mxCreateString( str );
	own = true;
	return *this;
}

object& object::operator=( const string &str )
{
	*this = str.data();
	return *this;
}

// Cast Operators -----------------------------------------------------------
//template<class T,class>
template<class T>
object::operator T() const
{
	mat<T> m = operator mat<T>();
	assert( m.size1() == 1 && m.size2() == 1 );
	return m(0,0);
}

template<class T>
object::operator vec<T>() const
{
	mat<T> m = operator mat<T>();
	assert( m.size1() == 1 || m.size2() == 1 );
	vec<T> v;
	if( m.size1() == 1 )
		v = m.row(0);
	else
		v = m.column(0);
	return v;
}

template<class T>
object::operator mat<T>() const
{
	cmat x   = operator cmat();
	mat<T> y = real(x);
	return y;
}

/// @cond INCLUDE_EXPLICIT_TEMPLATES
template<>
object::operator cmat() const
{
	// Get pointers
	mwSize ndim = mxGetNumberOfDimensions( (mxArray*)pdata );
	if( ndim != 2 )
		cerr << "object - wrong number of dimensions = " << ndim << endl;
	mwSize I      = mxGetDimensions( (mxArray*)pdata )[0];
	mwSize J      = mxGetDimensions( (mxArray*)pdata )[1];
	double *preal = mxGetPr( (mxArray*)pdata );
	double *pimag = NULL;
	if( mxIsComplex( (mxArray*)pdata ) )
		pimag = mxGetPi( (mxArray*)pdata );

	// Fill the array
	cmat m(I,J);
	for( uint i = 0; i < I; i++ )
		for( uint j = 0; j < J; j++ ) {
			vec<mwIndex> idx = {i,j};
			int offset = mxCalcSingleSubscript( (mxArray*)pdata, idx.size(), idx.data().begin() );
			m(i,j) = preal[offset];
			if( mxIsComplex( (mxArray*)pdata ) )
				m(i,j) += i*pimag[offset];
		}

	return m;
}
/// @endcond

void* object::data() const
{
	return pdata;
}

void object::data( void *p )
{
	if( own )
		mxDestroyArray( (mxArray*)pdata );
	pdata = p;
	own = true;
}

// Explicit instantiations -----------------------------------------------------------
/// @cond INCLUDE_EXPLICIT_TEMPLATES
template object& object::operator=( bool );
template object& object::operator=( int );
template object& object::operator=( double );
template object& object::operator=( const bvec& );
template object& object::operator=( const ivec& );
template object& object::operator=( const rvec& );
template object& object::operator=( const bmat& );
template object& object::operator=( const imat& );
template object& object::operator=( const rmat& );
template object::operator bool() const;
template object::operator int() const;
template object::operator double() const;
template object::operator complex() const;
template object::operator bvec() const;
template object::operator ivec() const;
template object::operator rvec() const;
template object::operator cvec() const;
template object::operator bmat() const;
template object::operator imat() const;
template object::operator rmat() const;
/// @endcond

//------------------------------------------------------------------------------
// Generic function call
//------------------------------------------------------------------------------
void trace( bool tf )
{
	_trace = tf;
}

}
}

#endif

/*
 * TI VOGLIO TANTO BENE :* (smak!)
 */
