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
// Expose
//------------------------------------------------------------

//int damfun( int a )
//{
//	return 2*a;
//}
//
//BOOST_PYTHON_MODULE( libbealab )
//{
//    using namespace boost::python;
//    def( "damfun", damfun );
//}

//------------------------------------------------------------
// Embed
//------------------------------------------------------------

static int py_init_flag = 0;

void _init()
{
	if(!py_init_flag){
		py_init_flag = 1;
		Py_Initialize();
		PyRun_SimpleString("import sys\n");
		PyRun_SimpleString("sys.path.append('.')\n");
		import_array();
	}
}

// Class members of object -----------------------------------
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
object::object( const Vec<T> &val )
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
template object::object( const Vec<bool>& );
template object::object( const Vec<int>& );
template object::object( const Vec<double>& );
template object::object( const Vec<complex>& );
/// @endcond

template<class T>
object::object( const Mat<T> &val )
{
	int N = val.size1();
	pdata = PyList_New( N );
	for( int n = 0; n < N; n++ ) {
		object o = Vec<T>( val.row( n ) );
		PyList_SetItem( (PyObject*)pdata, n, (PyObject*)o.data() );
		Py_INCREF( (PyObject*)o.data() );
	}
}
/// @cond INCLUDE_EXPLICIT_TEMPLATES
template object::object( const Mat<bool>& );
template object::object( const Mat<int>& );
template object::object( const Mat<double>& );
template object::object( const Mat<complex>& );
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
object::operator Vec<T>() const
{
	// Convert to PyArray
	PyObject* pd = PyArray_FROM_OTF( (PyObject*) pdata, NPY_DOUBLE, NPY_IN_ARRAY );
	if( pd == NULL )
		error("Cannot convert data to PyArray");
	double* data = (double*)PyArray_DATA( pd );
	if( data == NULL )
		error("Cannot get data pointer of PyArray");

	// Convert to Vec<T>
	int N = PyArray_DIM( pd, 0 );
	Vec<T> rv(N);
	for( int n = 0; n < N; n++ )
		rv(n) = data[n];
	return rv;
}
/// @cond INCLUDE_EXPLICIT_TEMPLATES
template object::operator Vec<int>() const;
template object::operator Vec<double>() const;
template object::operator Vec<complex>() const;
/// @endcond

template<class T>
object::operator Mat<T>() const
{
	// Convert to PyArray
	PyObject* pd = PyArray_FROM_OTF( (PyObject*) pdata, NPY_DOUBLE, NPY_IN_ARRAY );
	if( pd == NULL )
		error("Cannot convert data to PyArray");
	double* data = (double*)PyArray_DATA( pd );		// XXX Deberia ser complex
	if( data == NULL )
		error("Cannot get data pointer of PyArray");

	// Convert to Mat<T>
	int I = PyArray_DIM( pd, 0 );
	int J = PyArray_DIM( pd, 1 );
	Mat<T> rv(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			rv(i,j) = data[j+J*i];
	return rv;
}
/// @cond INCLUDE_EXPLICIT_TEMPLATES
template object::operator Mat<int>() const;
template object::operator Mat<double>() const;
template object::operator Mat<complex>() const;
/// @endcond

// Class members of functor ---------------------------------------------------

functor::functor( const string& module, const string& function )
{
	// Get function
    _init();
    PyObject *pName   = PyString_FromString( module.data() );
    PyObject *pModule = PyImport_Import( pName );
	if( pModule == NULL )
		error( "Cannot import module '" + module + "'" );
    pFunc   = PyObject_GetAttrString( pModule, function.data() );
	if( pFunc == NULL )
		error( "Module '" + module + "' has no member called '" + function + "'" );
}

object functor::_fcall( const deque<object> &O )
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

}
}

#endif

/*
 * TI VOGLIO TANTO BENE :* (smak!)
 */
