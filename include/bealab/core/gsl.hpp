// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/core/gsl.hpp
/// Interface to call GSL functions.

#ifndef _BEALAB_GSL_
#define	_BEALAB_GSL_

#include <bealab/core/blas/vector.hpp>
#include <bealab/core/blas/matrix.hpp>

namespace bealab
{
/// GSL interface
namespace gsl
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
