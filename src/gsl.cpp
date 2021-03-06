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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <bealab/core/blas.hpp>
#include <bealab/core/gsl.hpp>

namespace bealab
{
namespace gsl
{
//------------------------------------------------------------------------------
// Functor manager
//------------------------------------------------------------------------------
double sfunction_proxy( double x, void* pfun )
{
	function<double(double)> fun = *(function<double(double)>*)pfun;
	return fun( x );
}

double vfunction_proxy( double* px, size_t I, void* pfun )
{
	function<double(const rvec&)> fun = *(function<double(const rvec&)>*)pfun;
	rvec x(I);
	for( uint i = 0; i < I; i++ )
		x(i) = px[i];
	return fun( x );
}
//------------------------------------------------------------------------------
// GSL Vector interpreter class
//------------------------------------------------------------------------------
void vector::_alloc( int I )
{
	_data = gsl_vector_alloc( I );
	_own  = true;
}

void vector::_free()
{
	if(_own)
		gsl_vector_free( (gsl_vector*)_data );
	_own = false;
}

vector::vector( const rvec& x )
{
	int I = x.size();
	_alloc( I );
	*this = x;
}

vector& vector::operator=( const rvec& x )
{
	if( x.size() != ((gsl_vector*)_data)->size ) {
		cerr << x.size() << endl;
		cerr << ((gsl_vector*)_data)->size << endl;
		error( "gsl::vector: Sizes are different" );
	}
	int I = x.size();
	for ( int i = 0; i < I; i++ )
		gsl_vector_set( (gsl_vector*)_data, i, x(i) );
	return *this;
}

template<class T>
vector::operator vec<T>() const
{
	int I = ((gsl_vector*)_data)->size;
	vec<T> rv(I);
	for ( int i = 0; i < I; i++ )
		rv(i) = gsl_vector_get( (gsl_vector*)_data, i );
	return rv;
}

/// @cond INCLUDE_EXPLICIT_TEMPLATES
template vector::operator rvec() const;
/// @endcond

//------------------------------------------------------------------------------
// GSL Matrix interpreter class
//------------------------------------------------------------------------------
void matrix::_alloc( int I, int J )
{
	_data = gsl_matrix_alloc( I, J );
	_own  = true;
}

void matrix::_free()
{
	if(_own)
		gsl_matrix_free( (gsl_matrix*)_data );
	_own = false;
}

matrix::matrix( const rmat& x )
{
	int I = x.size1();
	int J = x.size2();
	_alloc( I, J );
	*this = x;
}

matrix& matrix::operator=( const rmat& x )
{
	if( x.size1() != ((gsl_matrix*)_data)->size1 || x.size2() != ((gsl_matrix*)_data)->size2 )
		error( "gsl::matrix: Sizes are different" );
	int I = x.size1();
	int J = x.size2();
	for ( int i = 0; i < I; i++ )
		for ( int j = 0; j < J; j++ )
			gsl_matrix_set( (gsl_matrix*)_data, i, j, x(i,j) );
	return *this;
}

template<class T>
matrix::operator mat<T>() const
{
	int I = ((gsl_matrix*)_data)->size1;
	int J = ((gsl_matrix*)_data)->size2;
	mat<T> rv(I,J);
	for ( int i = 0; i < I; i++ )
		for ( int j = 0; j < J; j++ )
			rv(i,j) = gsl_matrix_get( (gsl_matrix*)_data, i, j );
	return rv;
}

/// @cond INCLUDE_EXPLICIT_TEMPLATES
template matrix::operator rmat() const;
/// @endcond

}
}
