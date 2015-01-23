// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/core/blas/creation.hpp
/// Operations for creating vectors and matrices.

#ifndef _BEALAB_BLAS_CREATION_
#define	_BEALAB_BLAS_CREATION_

#include <bealab/core/blas/matrix.hpp>

namespace bealab
{
/// @defgroup blas_creation Creation
/// Operations for creating vectors and matrices.
/// @{


/// Construct a vector with the data given in a container of values
template<class value_type>
class vmake : public Vec<value_type> {

	template<class container>
	void vmake_( const container& cont )
	{
		this->resize( cont.size() );
		int i = 0;
		for( auto it = cont.begin(); it != cont.end(); it++ ) {
			(*this)(i) = *it;
			i++;
		}
	}

public:

	template<class E>
	vmake( const vector_interface<E>& v ) : Vec<value_type>(v) {}

	template<class container>
	vmake( const container& cont )
	{ vmake_(cont); }

	vmake( const initializer_list<value_type>& list )
	{ vmake_(list); }
};


/// Construct a vector by concatenating the sub-vectors in a container
template<class value_type>
class vconcat : public Vec<value_type> {

	template<class container>
	void vconcat_( const container& cont )
	{
		// Init
		auto p = cont.begin();
		int I  = cont.size();

		// Compute total size
		int L = 0;
	    for(int i = 0; i < I; i++)
	        L += p[i].size();

	    // Fill the vector
	    this->resize(L);
	    int cursor = 0;
	    for(int i = 0; i < I; i++) {
	    	int l = p[i].size();
	    	(*this)(range(cursor,cursor+l)) = p[i];
	    	cursor += l;
	    }
	}

public:

	template<class E>
	vconcat( const vector_interface<E>& cont )
	{ vconcat_<Vec<Vec<value_type>>>(cont); }

	vconcat( const initializer_list<vmake<value_type>> &list )
	{ vconcat_(list); }
};

/// Construct a matrix with the data given in nested initializer_list's
template<class value_type>
struct mmake : public Mat<value_type> {

	template<class E>
	mmake( const matrix_interface<E>& m ) : Mat<value_type>(m) {}

	mmake( const initializer_list<initializer_list<value_type>> &l )
	{
		// Compute sizes
		const initializer_list<value_type>* prow = l.begin();
		int	I = l.size();
		int	J = prow[0].size();
		for( int i = 1; i < I; i++ )
			assert( J == (int)prow[i].size() );

		// Fill the matrix
		this->resize(I,J);
		for( int i = 0; i < I; i++ )
			for( int j = 0; j < J; j++ )
				(*this)(i,j) = prow[i].begin()[j];
	}
};

/// Construct a matrix by concatenating the sub-matrices in nested
/// initializer_list's or a container
template<class value_type>
struct mconcat : public Mat<value_type> {

	mconcat( const initializer_list<initializer_list<mmake<value_type>>> &list )
	{
		// Compute sizes
		const initializer_list<mmake<value_type>>* prow = list.begin();
		int	Ib = list.size();
		int I  = 0;
		int J  = 0;
		for( int ib = 0; ib < Ib; ib++ ) {										// For each block row
			const mmake<value_type>* pblock = prow[ib].begin();					// Pointer to a block matrix in this row
			int Il = pblock[0].size1();											// Number of rows in this block row
			I     += Il;														// Update the total number of rows
			int Jb = prow[ib].size();											// Number of block columns in this row
			int Jl = 0;															// Total number of columns in this block row
			for( int jb = 0; jb < Jb; jb++ ) {
				assert( (int)pblock[jb].size1() == Il );
				Jl += pblock[jb].size2();
			}

			// Check that all block rows have the same number of columns
			if( ib == 0 )
				J = Jl;
			else
				assert( J == Jl );
		}

		// Fill the matrix
		this->resize(I,J);
		int cursor_I = 0;
		for( int ib = 0; ib < Ib; ib++ ) {
			int Jb = prow[ib].size();											// Number of block columns in this row
			const mmake<value_type>* pblock = prow[ib].begin();					// Pointer to a block matrix in this row
			int cursor_J = 0;
			int Il = pblock[0].size1();
			for( int jb = 0; jb < Jb; jb++ ) {
				int Jl = pblock[jb].size2();
				auto rI = range( cursor_I, cursor_I + Il );
				auto rJ = range( cursor_J, cursor_J + Jl );
				(*this)(rI,rJ) = pblock[jb];
				cursor_J += Jl;
			}
			cursor_I += Il;
		}
	}

	mconcat( const Mat<Mat<value_type>> &mmat )
	{
		// Compute sizes
		int	Ib = mmat.size1();
		int	Jb = mmat.size2();
		int I  = 0;
		int J  = 0;
		for( int ib = 0; ib < Ib; ib++ ) {										// For each block row
			int Il = mmat(ib,0).size1();										// Number of rows in this block row
			I     += Il;														// Update the total number of rows
			int Js = 0;															// Total number of columns in this block row
			for( int jb = 0; jb < Jb; jb++ ) {
				assert( (int)mmat(ib,jb).size1() == Il );
				Js += mmat(ib,jb).size2();
			}

			// Check that all block rows have the same number of columns
			if( ib == 0 )
				J = Js;
			else
				assert( J == Js );
		}

		// Fill the matrix
		this->resize(I,J);
		int cursor_I = 0;
		for( int ib = 0; ib < Ib; ib++ ) {
			int cursor_J = 0;
			int Il = mmat(ib,0).size1();
			for( int jb = 0; jb < Jb; jb++ ) {
				int Jl = mmat(ib,jb).size2();
				auto rI = range( cursor_I, cursor_I + Il );
				auto rJ = range( cursor_J, cursor_J + Jl );
				(*this)(rI,rJ) = mmat(ib,jb);
				cursor_J += Jl;
			}
			cursor_I += Il;
		}
	}
};


/// N-dimensional vector of zeros
inline
vector_interface<ublas::zero_vector<int>>
zeros( int N ) { return ublas::zero_vector<int>(N); }

/// N by M-dimensional matrix of zeros
inline
matrix_interface<ublas::zero_matrix<double>>
zeros( int N, int M ) {return ublas::zero_matrix<double>(N,M); }

/// N-dimensional vector of ones
inline
vector_interface<ublas::scalar_vector<int>>
ones( int N ) { return ublas::scalar_vector<int>(N); }

/// N by M-dimensional matrix of ones
inline
matrix_interface<ublas::scalar_matrix<double>>
ones( int N, int M ) {return ublas::scalar_matrix<double>(N,M); }

/// N-dimensional vector of zeros with a one in the i-th entry
inline
vector_interface<ublas::unit_vector<int>>
unit( int N, int i ) { return ublas::unit_vector<int>(N,i); }

/// N by M-dimensional identity matrix
inline
matrix_interface<ublas::identity_matrix<double>>
eye( int N ) {return ublas::identity_matrix<double>(N); }

/// Vector with values from 'from' to 'to_plus_1 - 1'
ivec vrange( int from, int to_plus_1 );

/// Vector with N values from 'from' by skipping 'skip' values
rvec vslice( double from, double skip, double N );

/// Vector with N values from 'from' to 'to', regularly spaced in the linear scale.
rvec linspace( double from, double to, int N );

/// Vector with N values from 'from' to 'to', regularly spaced in the logarithmic scale.
rvec logspace( double a, double b, int N );

/// @}
}
#endif
