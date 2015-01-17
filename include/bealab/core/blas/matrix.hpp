/// @file bealab/core/blas/matrix.hpp
/// Matrix modeling

#ifndef _BEALAB_BLAS_MATRIX_
#define	_BEALAB_BLAS_MATRIX_

#include <bealab/core/blas/vector.hpp>

namespace bealab
{
/// @defgroup blas_matrix Matrices
/// Matrix modeling
/// @{

/// Matrix class with extended construction
template<class value_type>
//class matrixx : public ublas::matrix<value_type,ublas::column_major> {
class matrixx : public ublas::matrix<value_type> {
public:

	using ublas::matrix<value_type>::matrix;
	using ublas::matrix<value_type>::operator=;

	/// Default constructor
	matrixx() = default;

	/// Construct a matrix with the data given in two nested
	/// initializer_list's
	matrixx( const initializer_list<initializer_list<value_type>> &l )
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

	/// Construct a matrix by concatenating the sub-matrices given in two
	/// nested initializer_list's
	matrixx( const initializer_list<initializer_list<matrixx<value_type>>> &list )
	{
		// Compute sizes
		const initializer_list<matrixx<value_type>>* prow = list.begin();
		int	Ib = list.size();
		int I  = 0;
		int J  = 0;
		for( int ib = 0; ib < Ib; ib++ ) {										// For each block row
			const matrixx<value_type>* pblock = prow[ib].begin();				// Pointer to a block matrix in this row
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
			const matrixx<value_type>* pblock = prow[ib].begin();				// Pointer to a block matrix in this row
			int cursor_J = 0;
			int Il = pblock[0].size1();
			for( int jb = 0; jb < Jb; jb++ ) {
				int Jl = pblock[jb].size2();
				auto rI = range( cursor_I, cursor_I + Il );
				auto rJ = range( cursor_J, cursor_J + Jl );
				ublas::matrix_range<matrixx<value_type>>(*this,rI,rJ)
						= pblock[jb];
				cursor_J += Jl;
			}
			cursor_I += Il;
		}
	}

	/// Copy into the matrix the data given in two nested initializer_list's
	matrixx<value_type>& operator=(
			const initializer_list<initializer_list<value_type>> &list )
	{
		*this = matrixx<value_type>(list);
		return *this;
	}

	/// Copy into the matrix the data obtained by concatenating the sub-matrices
	/// given in two nested initializer_list's
	matrixx<value_type>& operator=(
			const initializer_list<initializer_list<matrixx<value_type>>> &list )
	{
		*this = matrixx<value_type>(list);
		return *this;
	}
};

/// Interface for any matrix_expression
template<class base>
class matrix_interface : public base {

	/// Compute a range ublas vector to implement row/column indirect access
	static
	ublas::vector<int> _vrange( int a, int b )
	{
		int N = b - a;
		ublas::vector<int> v(N);
		for( int n = 0; n < N; n++ )
			v(n) = a+n;
		return v;
	}

public:

	using typename base::value_type;
	using typename base::const_reference;
	using typename base::reference;
	using base::base;
	using base::operator=;

	/// @name Constructors

	/// Default constructor (because it is not inherited)
	matrix_interface() = default;

	/// Copy constructor using the base
	matrix_interface( const base& b ) : base(b) {}
	/// @}

	/// @name Base access

	/// Access the base (constant)
	const base& operator()() const
	{
		return *this;
	}

	/// Access the base (non-constant)
	base& operator()()
	{
		return *this;
	}
	///@}

	/// @name Entry-wise access

	/// Access an entry (constant)
	const_reference operator()( int a, int b ) const
	{
		return base::operator()( a, b );
	}

	/// Access an entry (non-constant)
	reference operator()( int a, int b )
	{
		return base::operator()( a, b );
	}

	/// Access a two-dimensional range
	matrix_interface<ublas::matrix_range<base>>
	operator()( const range& r1, const range& r2 ) const
	{
		auto& self = const_cast<matrix_interface<base>&>(*this);
		return { self, r1, r2 };
	}

	/// Access a two-dimensional slice
	matrix_interface<ublas::matrix_slice<base>>
	operator()( const slice& s1, const slice& s2 ) const
	{
		auto& self = const_cast<matrix_interface<base>&>(*this);
		return { self, s1, s2 };
	}

	/// Indirect two-dimensional access
	matrix_interface<ublas::matrix_indirect<base>>
	operator()( const indirect& i1,
				const indirect& i2 ) const
	{
		auto& self = const_cast<matrix_interface<base>&>(*this);
		return { self, i1, i2 };
	}
	/// @}

	/// @name Row-wise access

	/// Access a row
	vector_interface<ublas::matrix_row<base>> row( int r ) const
	{
		auto& self = const_cast<matrix_interface<base>&>(*this);
		return { self, (unsigned int)r };
	}

	/// Access a row range
	matrix_interface<ublas::matrix_range<base>> row( const range& r1 ) const
	{
		auto& self = const_cast<matrix_interface<base>&>(*this);
		return { self, r1, range(0,this->size2()) };
	}

	/// Access a row slice
	matrix_interface<ublas::matrix_slice<base>> row( const slice& s1 ) const
	{
		auto& self = const_cast<matrix_interface<base>&>(*this);
		return { self, s1, slice(0,1,this->size2()) };
	}

	/// Indirect row access
	matrix_interface<ublas::matrix_indirect<base>> row( const indirect& i1 ) const
	{
		auto& self = const_cast<matrix_interface<base>&>(*this);
		return { self, i1, indirect(_vrange(0,this->size2())) };
	}
	/// @}

	/// @name Column-wise access

	/// Access a column
	vector_interface<ublas::matrix_column<base>> column( int c ) const
	{
		auto& self = const_cast<matrix_interface<base>&>(*this);
		return { self, (unsigned int)c };
	}

	/// Access a column range
	matrix_interface<ublas::matrix_range<base>> column( const range& r2 ) const
	{
		auto& self = const_cast<matrix_interface<base>&>(*this);
		return { self, range(0,this->size1()), r2 };
	}

	/// Access a column slice
	matrix_interface<ublas::matrix_slice<base>> column( const slice& s2 ) const
	{
		auto& self = const_cast<matrix_interface<base>&>(*this);
		return { self, slice(0,1,this->size1()), s2 };
	}

	/// Indirect column access
	matrix_interface<ublas::matrix_indirect<base>> column( const indirect& i2 ) const
	{
		auto& self = const_cast<matrix_interface<base>&>(*this);
		return { self, indirect(_vrange(0,this->size1())), i2 };
	}
	/// @}
};

/// Dense matrix template
template<class value_type>
using Mat = matrix_interface<matrixx<value_type>>;

/// Triangular matrix template
template<class value_type, class LU>
using triangular_matrix = matrix_interface<
		ublas::triangular_matrix<value_type,LU>>;

/// Symmetric matrix template
template<class value_type>
using symmetric_matrix = matrix_interface<ublas::symmetric_matrix<value_type>>;

/// Hermitian matrix template
template<class value_type>
using hermitial_matrix = matrix_interface<ublas::hermitian_matrix<value_type>>;

/// Banded matrix template
template<class value_type>
using banded_matrix = matrix_interface<ublas::banded_matrix<value_type>>;

/// Sparse matrix template
template<class value_type>
using sparse_matrix = matrix_interface<ublas::compressed_matrix<value_type>>;

// Shorthand expressions for dense matrices
typedef	Mat<bool> bmat;															///< Boolean dense matrix
typedef	Mat<int> imat;															///< Integer dense matrix
typedef Mat<double> rmat;														///< Real dense matrix
typedef	Mat<complex> cmat;														///< Complex dense matrix

/// Convert an expression template into a temporary
template<class E>
matrix_interface<matrixx<typename E::value_type>>
noproxy( const matrix_interface<E>& x )
{
	return x;
}

/// Displays a matrix in the console
template<class A, class B, class T>
std::basic_ostream<A,B> &operator<<( std::basic_ostream<A,B> &os,
                                     const matrix_interface<T> &m )
{
    int I = m.size1();
    int J = m.size2();
	os << '[';
    for( int i = 0; i < I; i++ ) {
    	os << '[';
    	if( J > 0 )
        	os << m(i,0);
        for( int j = 1; j < J; j++ )
        	os << ", " << m(i,j);
    	os << ']';
    	if( i != I-1 )
    		os << endl;
    }
	os << ']';
    return os;
}

/// Apply a function to all the elements of a matrix
template<class F, class T>
Mat<typename result_of<F(typename T::value_type)>::type>
	entrywise( F fun, const matrix_interface<T>& A )
{
	typedef typename result_of<F(typename T::value_type)>::type R;
	int I = A.size1();
	int J = A.size2();
	Mat<R> B(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			B(i,j) = fun( A(i,j) );
	return B;
}

/// @}
}
#endif
