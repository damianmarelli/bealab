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
#ifndef BEALAB_NOPYTHON

#include <python2.7/Python.h>
#include <python2.7/numpy/arrayobject.h>
#include <bealab/core/python.hpp>

namespace bealab
{
namespace python
{
//------------------------------------------------------------
// Object
//------------------------------------------------------------
object::object( const object &o )
{
	pdata = o.pdata;
	Py_INCREF( pdata );
}

object::object( int val )
{
	pdata = PyInt_FromLong( val );
}

object::object( double val )
{
	pdata = PyFloat_FromDouble( val );
}

object::object( complex val )
{
	pdata = PyComplex_FromDoubles( val.real(), val.imag() );
}

object::object( const char *val )
{
	pdata = PyString_FromString( val );
}

object::object( const string &val )
{
	pdata = PyString_FromString( val.data() );
}

template<class T>
object::object( const vec<T> &val )
{
	int N = val.size();
	pdata = PyList_New( N );
	for( int n = 0; n < N; n++ ) {
		object o = val(n);
		PyList_SetItem( (PyObject*)pdata, n, (PyObject*)o.data() );
		Py_INCREF( (PyObject*)o .data());
	}
}
/// @cond INCLUDE_EXPLICIT_TEMPLATES
template object::object( const vec<bool>& );
template object::object( const vec<int>& );
template object::object( const vec<double>& );
template object::object( const vec<complex>& );
/// @endcond

template<class T>
object::object( const mat<T> &val )
{
	int N = val.size1();
	pdata = PyList_New( N );
	for( int n = 0; n < N; n++ ) {
		object o = vec<T>( val.row( n ) );
		PyList_SetItem( (PyObject*)pdata, n, (PyObject*)o.data() );
		Py_INCREF( (PyObject*)o.data() );
	}
}
/// @cond INCLUDE_EXPLICIT_TEMPLATES
template object::object( const mat<bool>& );
template object::object( const mat<int>& );
template object::object( const mat<double>& );
template object::object( const mat<complex>& );
/// @endcond

object::~object()
{
	if( pdata != NULL ) {
//			cout << pdata->ob_refcnt << endl;
		Py_DECREF(pdata);
	}
}

void object::operator=( const object &o )
{
	// Destruct
	this->~object();

	// Reconstruct
	pdata = o.pdata;
	Py_INCREF( pdata );
}

object::operator int() const
{
	long rv = PyInt_AsLong( (PyObject*)pdata );
	if( PyErr_Occurred() )
		error("Error converting Python object to int");
	return rv;
}

object::operator double() const
{
	double rv = PyFloat_AsDouble( (PyObject*)pdata );
	if( PyErr_Occurred() )
		error("Error converting Python object to double");
	return rv;
}

object::operator complex() const
{
	double real = PyComplex_RealAsDouble( (PyObject*)pdata );
	if( PyErr_Occurred() )
		error("Error converting Python object to complex");
	double imag = PyComplex_ImagAsDouble( (PyObject*)pdata );
	if( PyErr_Occurred() )
		error("Error converting Python object to complex");
	return complex(real,imag);
}

object::operator string() const
{
	string rv = PyString_AsString( (PyObject*)pdata );
	if( PyErr_Occurred() )
		error("Error converting Python object to string");
	return rv;
}

template<class T>
object::operator vec<T>() const
{
	// Convert to PyArray
	PyObject* pd = PyArray_FROM_OTF( (PyObject*) pdata, NPY_DOUBLE, NPY_IN_ARRAY );
	if( pd == NULL )
		error("Cannot convert data to PyArray");
	double* data = (double*)PyArray_DATA( pd );
	if( data == NULL )
		error("Cannot get data pointer of PyArray");

	// Convert to vec<T>
	int N = PyArray_DIM( pd, 0 );
	vec<T> rv(N);
	for( int n = 0; n < N; n++ )
		rv(n) = data[n];
	return rv;
}
/// @cond INCLUDE_EXPLICIT_TEMPLATES
template object::operator vec<int>() const;
template object::operator vec<double>() const;
template object::operator vec<complex>() const;
/// @endcond

template<class T>
object::operator mat<T>() const
{
	// Convert to PyArray
	PyObject* pd = PyArray_FROM_OTF( (PyObject*) pdata, NPY_DOUBLE, NPY_IN_ARRAY );
	if( pd == NULL )
		error("Cannot convert data to PyArray");
	double* data = (double*)PyArray_DATA( pd );		// XXX Deberia ser complex
	if( data == NULL )
		error("Cannot get data pointer of PyArray");

	// Convert to mat<T>
	int I = PyArray_DIM( pd, 0 );
	int J = PyArray_DIM( pd, 1 );
	mat<T> rv(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			rv(i,j) = data[j+J*i];
	return rv;
}
/// @cond INCLUDE_EXPLICIT_TEMPLATES
template object::operator mat<int>() const;
template object::operator mat<double>() const;
template object::operator mat<complex>() const;
/// @endcond

//------------------------------------------------------------
// Function
//------------------------------------------------------------

int _function_handler::py_init_flag = 0;

void _function_handler::init()
{
	if(!py_init_flag){
		py_init_flag = 1;
		Py_Initialize();
		PyRun_SimpleString("import sys\n");
		PyRun_SimpleString("sys.path.append('.')\n");
		import_array();
	}
}

object _function_handler::fcall( const deque<object> &O )
{
	// Make arguments
	int N = O.size();
	PyObject *pArgs = PyTuple_New( N );
	for( int i = 0; i < N; ++i) {
		PyTuple_SetItem( pArgs, i, (PyObject*)O[i].data() );
		Py_INCREF( (PyObject*)O[i].data() );
	}

	// Call function
	object rv = PyObject_CallObject( (PyObject*)pFunc, pArgs );
	if( PyErr_Occurred() ) {
		PyErr_Print();
		error("");
	}

	// Finalize
	Py_DECREF(pArgs);

    return rv;
}

_function_handler::_function_handler( const string& module, const string& function )
{
	// Get function
    init();
    PyObject *pName   = PyString_FromString( module.data() );
    PyObject *pModule = PyImport_Import( pName );
	if( pModule == NULL )
		error( "Cannot import module '" + module + "'" );
    pFunc   = PyObject_GetAttrString( pModule, function.data() );
	if( pFunc == NULL )
		error( "Module '" + module + "' has no member called '" + function + "'" );
}

}
}

#endif

/*
 * TI VOGLIO TANTO BENE :* (smak!)
 */
