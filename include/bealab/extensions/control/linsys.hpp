// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/extensions/control/linsys.hpp
/// Linear system models.

#ifndef _BEALAB_LINSYS_
#define	_BEALAB_LINSYS_

#include <bealab/scilib/sequences.hpp>
#include <bealab/scilib/linalg.hpp>
#include <bealab/scilib/stats.hpp>

namespace bealab
{
namespace control
{
//------------------------------------------------------------------------------
/// @defgroup linsys Linear systems
/// Linear system models.
/// @{

/// ARMA model.
/// The template parameter S is the type of input/output and T is the type of
/// each entry of the numerator and denominator vectors.
template<class T, class S=T>
class arma {

//	struct _ST {
//		vec<S> input;
//		vec<S> output;
//		int ibuflength, obuflength;
////		_ST() : ibuflength(0), ibuflength(0) {}
//		void shift( vec<S>& buffer, const S& x )
//		{
//			int N = buffer.size();
//			for( int n = N-1; n > 0; n-- )
//				buffer(n) = buffer(n-1);
//			if( N > 0 )
//			buffer(0) = x;
//		}
//	};

	/// Internal state structure
	class buffer_t : public vec<S> {
		int buflen;
	public:
		buffer_t() : vec<S>(), buflen(0) {}

		void set_len( int len )
		{
			if( len < buflen )
				this->resize(len);
			buflen = len;
		}

		void push( const S& x )
		{
			if( (int)this->size() < buflen )
				this->resize(this->size()+1);
			int N = this->size();
			for( int n = N-1; n > 0; n-- )
				(*this)(n) = (*this)(n-1);
			if( N > 0 )
				(*this)(0) = x;
		}

		S average( const vec<T>& x )
		{
			int N = this->size();
			const vec<S>& prod = element_prod( x(range(0,N)), *this );
			return sum( prod );
		}
	};

	/// Internal state
//	_ST _state;
	buffer_t inbuf, outbuf;

	vec<T> _b;	///< Numerator
	vec<T> _a0;	///< First element of the denominator. Must be 1.
	vec<T> _ar;	///< Denominator remainder (= denominator - 1)

	/// Return type
//	typedef decltype(T()*S()) R;

public:

	/// @name Constructors
	arma() { clear(); }
	arma( const vec<T> &b_, const vec<T> &a_ )
	{
		set_coeffs( b_, a_ );
	}
	arma( const arma<T,S> &filt ) { *this = filt; }
	/// @}

	/// @name Get and set

	/// Clear state
	void clear()
	{
		inbuf.clear();
		outbuf.clear();
	}

	/// Cppy operator
	void operator=( const arma<T,S> &filt )
	{
		set_coeffs( filt.num(), filt.den() );
		clear();
	}

	/// Get numerator
	vec<T> num() const { return _b; }

	/// Get Denominator
	vec<T> den() const { return {_a0, _ar}; }

	/// Set numerator and denominator
	void set_coeffs( const vec<T> &b_, const vec<T> &a_ )
	{
		// Set b
		_b  = b_;

		// Set a
		if( a_.size() == 0 ) {
			_a0.resize(0);
			_ar.resize(0);
		}
		else {
			if( norm( a_(0)*a_(0) - a_(0) ) != 0 )
				error("arma<T,S>.set_coeffs() - The first denominator element is not normalized");
			_a0 = a_( range(0,1) );
			_ar = a_( range( 1, a_.size() ) );
		}

		// Set buffers
		inbuf.set_len( _b.size() );
		outbuf.set_len( _ar.size() );
	}
	/// @}

	/// @name Filtering

	/// Filter one sample
	S operator()( const S &x )
	{
		inbuf.push( x );
//		R y = sum(element_prod( _b, _state.input ))
//		    - sum(element_prod( _ar, _state.output ));
//		const vec<S>& Y1 = element_prod( _b, _state.input );
//		const vec<S>& Y2 = element_prod( _ar, _state.output );
//		const S& y1 = sum(Y1);
//		const S& y2 = sum(Y2);
		S y  = inbuf.average( _b );
		if( outbuf.size() > 0 )
			y -= outbuf.average( _ar );
		outbuf.push( y );
		return y;
	}

	/// Filter a sequence
	sequence<S> operator()( const sequence<S> &x )
	{
		sequence<S> y( x.size(), x.t1() );
		for( int t = x.t1(); t <= x.t2(); t++ )
			y(t) = (*this)( x(t) );
		return y;
	}
	/// @}

	/// @name Response

	/// Impulse response
	sequence<S> impulse_response( int N ) const
	{
		arma<T,S> local = *this;
		vec<S> x_ = { {1}, zeros(N-1) };
		sequence<S> x( x_, 0  );
		return local( x );
	}

	/// Frequency response
	cvec frequency_response( const rvec &W ) const
	{
		int I = W.size();
		cvec N(I), D(I);
		vec<T> a = {_a0,_ar};
		for( int i = 0; i < I; i++ ) {
			N(i) = polyval( _b, exp(j*W(i)) );
			D(i) = polyval( a, exp(j*W(i)) );
		}
		return element_div(N,D);
	}
	/// @}
};

/// Transfer function.
typedef arma<double> transfer_function;

/// State space model.
class state_space {

	bool noise_f = false;
	rmat Qh, Rh;

public:

	rmat A, B, C, D;
	rmat Q, R;
	rvec x;
	state_space() = default;
	state_space( const rmat &A_, const rmat &B_, const rmat &C_, const rmat& D_	) :
		A(A_), B(B_), C(C_), D(D_)
	{
		assert( A.size1() == A.size2() );
		assert( A.size2() == B.size1() );
		assert( C.size2() == A.size1() );
		assert( B.size2() == D.size2() );
		assert( C.size1() == D.size1() );
		reset();
	}

	state_space( const rmat &A_, const rmat &B_, const rmat &C_, const rmat& D_,
			const rmat& Q_, const rmat& R_ ) :
				state_space(A_,B_,C_,D_)
	{
		assert( Q_.size1() == Q_.size2() );
		assert( Q_.size1() == A.size1() );
		assert( R_.size1() == R_.size2() );
		assert( R_.size1() == C.size1() );
		Q  = Q_;
		R  = R_;
		Qh = real(msqrt(Q));
		Rh = real(msqrt(R));
		noise_f = true;
	}

	void reset()
	{
		x = zeros(A.size1());
	}

	rvec operator()( const rvec& u )
	{
		rvec y = C * x + D * u;
		x = A * x + B * u;
		if( noise_f ) {
			x += Qh * randn(Qh.size2());
			y += Rh * randn(Rh.size2());
		}
		return y;
	}

	sequence<rvec> operator()( const sequence<rvec>& u )
	{
		int T = u.size();
		vec<rvec> yvec(T);
		for( int t = 0; t < T; t++ )
			yvec(t) = (*this)( u.buffer()(t) );
		return { yvec, u.t1() };
	}

	rseq operator()( const rseq& u )
	{
		int T = u.size();
		rvec yvec(T);
		for( int t = 0; t < T; t++ )
			yvec(t) = (*this)( rvec{u.buffer()(t)} )(0);
		return { yvec, u.t1() };
	}

	sequence<rvec> impulse_response( int T )
	{
		vec<rvec> yvec(T);
		int N   = x.size();
		yvec(0) = (*this)( (rvec)ones(N) );
		for( int t = 1; t < T; t++ )
			yvec(t) = (*this)( (rvec)zeros(N) );
		return { yvec, 0 };
	}
};

state_space discretize_sampling( const state_space& ss, double T );

state_space controllable_canonical_form( const transfer_function& tf );

/// @name Support for arrays of transfer functions

/// Converts a matrix transfer function into a matrix of scalar transfer functions
mat<transfer_function> tfm2mtf( const arma<rmat,rvec> &S );

/// Converts a matrix of scalar transfer functions into a matrix transfer function
arma<rmat,rvec> mtf2tfm( const mat<transfer_function> &M );
/// @}

/// @}
}
}
#endif
