// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/scilib/arrays.hpp
/// Some functions to deal with arrays of objects.

#ifndef _BEALAB_ARRAYS_
#define	_BEALAB_ARRAYS_

#include <bealab/scilib/fourier.hpp>

namespace bealab
{

/// @defgroup arrays Arrays
/// Some functions to deal with arrays of objects.
/// @{

/// @name Conversions Mat/Vec
template<class T>
Mat<Vec<T>> vm2mv( const Vec<Mat<T>>& x )
{
	if( x.size() == 0 )
		return Mat<Vec<T>>();
	Mat<T> x0 = x(0);
	int I     = x0.size1();
	int J     = x0.size2();
	int N     = x.size();
	Mat<Vec<T>> y(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ ) {
			y(i,j).resize(N);
			for( int n = 0; n < N; n++ )
				y(i,j)(n) = x(n)(i,j);
		}
	return y;
}

template<class T>
Vec<Mat<T>> mv2vm( const Mat<Vec<T>>& x )
{
	int I = x.size1();
	int J = x.size2();
	if( x.size1() * x.size2() == 0 )
		return Mat<Vec<T>>(I,J);
	int N = x(0,0).zeros(I,J);
	Vec<Mat<T>> y(N);
	for( int n = 0; n < N; n++ ) {
		y(n).resize(I,J);
		for( int i = 0; i < I; i++ )
			for( int j = 0; j < J; j++ )
					y(n)(i,j) = x(i,j)(n);
	}
	return y;
}
/// @}

/// @name Conversions Vec/Seq
template<class T>
Vec<Seq<T>> sv2vs( const Seq<Vec<T>> &x )
{
	Vec<T> x1 = x(x.t1());
	int N = x1.size();
	Vec<Seq<T>> y(N);
	for( int i = 0; i < N; i++ ) {
		auto component = [i](const Vec<T>& x){return x(i);};
		y(i) = entrywise(component)( x );
	}
	return y;
}

template<class T>
Seq<Vec<T>> vs2sv( const Vec<Seq<T>> &x )
{
	int N = x.size();
	ivec T1 = zeros(N);
	ivec T2 = zeros(N);
	for( int i = 0; i < N; i++ ) {
		T1(i) = x(i).t1();
		T2(i) = x(i).t2();
	}
	int t1 = min(T1);
	int t2 = max(T2);

	Seq<Vec<T>> y( t2-t1+1, t1 );
	for( int t = t1; t <= t2; t++ ) {
		auto sample = [t](const Seq<T>& xi){return xi(t);};
		y(t) = entrywise(sample)(x);
	}
	return y;
}
/// @}

/// @name Conversions Mat/Seq
template<class T>
Mat<Seq<T>> sm2ms( const Seq<Mat<T>> &x )
{
	Mat<T> x1 = x(x.t1());
	int I = x1.size1();
	int J = x1.size2();
	Mat<Seq<T>> y(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ ) {
			auto component = [i,j](const Mat<T>& x){return x(i,j);};
			y(i,j) = entrywise(component)( x );
		}
	return y;
}

template<class T>
Seq<Mat<T>> ms2sm( const Mat<Seq<T>> &x )
{
	int I = x.size1();
	int J = x.size2();
	imat T1 = zeros(I,J);
	imat T2 = zeros(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ ) {
			T1(i,j) = x(i,j).t1();
			T2(i,j) = x(i,j).t2();
		}
	int t1 = min(T1);
	int t2 = max(T2);

	Seq<Mat<T>> y( t2-t1+1, t1 );
	for( int t = t1; t <= t2; t++ ) {
		auto sample_t = [t](const Seq<T>& xij){return xij(t);};
		y(t) = entrywise(sample_t)( x );
	}
	return y;
}
/// @}

/// @name Sampling
template<class T>
Vec<T> sample_get( const Vec<Seq<T> > &v, int t )
{
	int	I = v.size();
	Vec<T>	w(I);
	for( int i = 0; i < I; i++ )
		w(i) = v(i)(t);
	return w;
}

template<class T>
void sample_set( Vec<Seq<T> > &v, int t, const Vec<T> &w )
{
	int	I = v.size();
	for( int i = 0; i < I; i++ )
		v(i)(t) = w(i);
}

template<class T>
Mat<T> sample_get( const Mat<Seq<T>> &m, int t )
{
	int	I = m.size1();
	int	J = m.size2();
	Mat<T>	M(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			M(i,j) = m(i,j)(t);
	return M;
}

template<class T>
void sample_set( Mat<Seq<T> > &m, int t, const Mat<T> &M )
{
	int	I = m.rows();
	int	J = m.cols();
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			m(i,j)(t) = M(i,j);
}
/// @}

/// @name Inversion
template<class T>
Mat<Seq<T>> inv( const Mat<Seq<T>>& G, double tol=1e-12 )
{
	int N = max( G.apply( [](const cseq& x){ return x.size(); } ) );			// Maximum sequence size in G
	Mat<Seq<T>> Gi;
	for(;; N *= 2 ) {

		// Compute the inverse
		Mat<cvec> Gf   = G.apply( [N](const cseq& x){ return dtft(x,N); } );
		Vec<cmat> Gf_  = ms2sm( Mat<cseq>(Gf) );
		Vec<cmat> Gfi_ = element_inv( Gf_ );
		Mat<cvec> Gfi  = sm2ms( Seq<cmat>(Gfi_) );
		Gi             = apply( [N](const cvec& x){ return idtft(x); } )( Gfi );

		// Evaluate
		double err = norm( G - G*Gi*G ) / norm( G );
		cout << "inv(Mat<cseq>) - error = " << err << endl;
		if( err < tol )
			break;
	}

	return Gi;
}
/// @}

/// @}
}
#endif
