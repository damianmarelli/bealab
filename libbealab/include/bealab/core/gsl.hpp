/// @file bealab/core/gsl.hpp
/// Interface to call GSL functions.

#ifndef _BEALAB_GSL_
#define	_BEALAB_GSL_

#include <bealab/core/blas/vector.hpp>
#include <bealab/core/blas/matrix.hpp>

namespace bealab
{
namespace _gsl
{

/// @defgroup interfaces_gsl GSL interface
/// Interface to call GSL functions.
/// @{

/// Permits passing a scalar std::function to a GSL function.
/// To pass the functor fun to a GSL function do: gsl_function( functor_manager, &fun );
double sfunction_proxy( double x, void* pfun );

/// Permits passing a vector std::function to a GSL function.
double vfunction_proxy( double* px, size_t I, void* pfun );

/// Interpreter class for GSL vectors
class vector {

	void *_data;
	bool _own;
	void _alloc( int I );
	void _free();

public:
	/// @name Constructors & destructor
//	vector() : _data(NULL), _own(false) {}
	vector( int I ) : _own(true) { _alloc(I); }
	vector( const rvec& x );
	vector( const void *pgslvec ) : _data(const_cast<void*>(pgslvec)), _own(false) {}
	~vector() { _free(); }
	/// @}

	/// @name Assignment
	vector& operator=( const rvec& x );
	/// @}

	/// @name Cast operators
	template<class T>
	operator Vec<T>() const;
	template<class T>
	operator T*() const { return (T*)_data; }
	/// @}
};

/// Interpreter class for GSL matrices
class matrix {

	void *_data;
	bool _own;
	void _alloc( int I, int J );
	void _free();

public:
	/// @name Constructors & destructor
	matrix( int I, int J ) : _own(true) { _alloc(I,J); }
	matrix( const rmat& x );
	matrix( const void *pgslmat ) : _data(const_cast<void*>(pgslmat)), _own(false) {}
	~matrix() { _free(); }
	/// @}

	/// @name Assignment
	matrix& operator=( const rmat& x );
	/// @}

	/// @name Cast operators
	template<class T>
	operator Mat<T>() const;
	template<class T>
	operator T*() const { return (T*)_data; }
	/// @}
};

/// @}
}
}
#endif
