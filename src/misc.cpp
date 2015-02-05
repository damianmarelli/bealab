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
#include <bealab/core/python.hpp>
#include <bealab/core/matlab.hpp>
#include <bealab/scilib/fourier.hpp>
#include <bealab/scilib/arrays.hpp>
#include <bealab/extensions/signal/misc.hpp>

namespace bealab
{
namespace signal
{

#ifndef BEALAB_NOMATLAB
//state_space spectral_realization( const mat<rseq> &Rx )
//{
//	// 3D array for matlab call
//	sequence<rmat> Rx1 = ms2sm( Rx );
//	Rx1   = Rx1.trunc( 0, Rx1.t2() );
//	int I = Rx.size1();
//	int J = Rx.size2();
//	int K = Rx1.size();
//	ivec dims = { I, J, K };
//	double a[I][J][K];
//	for( int i = 0; i < I; i++ )
//		for( int j = 0; j < J; j++ )
//			for( int k = 0; k < K; k++ )
//				a[i][j][k] = Rx1(k)(i,j);
//
//	// Matlab call
//	matlab::object b( (double*)a, dims );
//	auto rv = matlab::functor<3>( "spectral_realization" )( b, 1e-3 );
//	rmat A = rv[0];
//	rmat B = rv[1];
//	rmat C = rv[2];
//
//	// Return ss model
//	return state_space( A, B, C, zeros(C.size1(),B.size2()) );
//}
#endif

mat<cseq> spectral_factorization( const mat<cseq>& X, double tol )
{
	class {
		int rankred( const rvec& x, double tol )
		{
			int N = x.size();
			for( int n = 0; n < N; n++ )
				if( x(n) < tol )
					return n;
			return N;
		}
	public:
		mat<cseq> operator()( const mat<cseq>& X, double tol, int N )
		{
			const mat<cvec>& Xf  = entrywise( [N](const cseq& s) {return dtft(s,N);} )( X );
			const vec<cmat>& Xf_ = ms2sm( mat<cseq>(Xf) ).buffer();
			vec<cmat> Yf_(N);
			vec<cmat> U(N);
			vec<rvec> d(N);
			double max_sv = 0;
			for( int n = 0; n < N; n++ ) {
				auto UDV = svd( Xf_(n) );
				U(n)   = get<0>( UDV );
				rmat D = get<1>( UDV );
				cmat V = get<2>( UDV );
				d(n) = sqrt(diag(D));
				max_sv = max( max_sv, max(d(n)) );
			}
			ivec r(N);
			for( int n = 0; n < N; n++ )
				r(n) = rankred( d(n), tol*max_sv );
			int R = max(r);
			for( int n = 0; n < N; n++ ) {
//				Yf_(n) = U(n).col_range( 0, R-1 ) * diag( d(n)(range(0,R)) );
				cmat tmp = U(n)( range(0,U(n).size1()), range(0,R) );
				Yf_(n)   = tmp * diag( d(n)(range(0,R)) );
			}
			const mat<cseq>& Yf = sm2ms( sequence<cmat>(Yf_) );
			const mat<cseq>& Y  = entrywise( [](const cseq& v) {return idtft(v.buffer());} )( Yf );
			return Y;
		}
	} spectral_factorization_N;

	const cseq& x = sum(X);
	int N         = x.size();
	mat<cseq> Y;
	double err    = inf;
	for( int k = 0; err > tol; k++ ) {
		Y   = spectral_factorization_N( X, tol/10, N*pow(2,k) );
		err = norm( X - Y * adjoint(Y) );
//		cout << "error = " << err << endl;
	}
	return Y;
}

//------------------------------------------------------------------------------
// FIR Filter design
//------------------------------------------------------------------------------
#ifndef BEALAB_NOPYTHON
control::transfer_function butter( int order, double bandwidth, bool analog )
{
	python::function<rmat> butter( "scipy.signal", "butter" );
	rmat r  = butter( order, bandwidth, "low", analog );
	rvec b  = r.row(0);
	rvec a  = r.row(1);
	return control::transfer_function( b, a );
}
#endif

double raisedcosine( double t, double Fs, double beta )
{
	double T = 1 / Fs;
	if( abs(2*beta/T * t) == 1 )
		return sinc( t/T ) * pi / 4;
	else
		return sinc( t/T ) * cos( pi * beta  / T * t ) /
			( 1 - pow( 2*beta/T * t, 2 ) );
}

rseq raisedcosine( int N, double ws, double beta )
{
	assert( abs(ws) <= pi );
	int n1    = -floor(N/2.);
	double Fs = ws/(2*pi);
	rvec v(N);
	for( int n = 0; n < N; n++ ) {
		double t = n + n1;
		v(n)     = raisedcosine( t, Fs, beta );
	}
	return { Fs*v, n1 };
}

//rseq raisedcosine( int N, double ws, double beta )
//{
//	int t1   = -floor(N/2.);
//	double T = 2*pi / ws;
//	rvec v(N);
//	for( int i = 0; i < N; i++ ) {
//		double t = i + t1;
//		if( abs(2*beta/T * t) == 1 )
//			v(i) = sinc( t/T ) * pi / 4;
//		else
//			v(i) = sinc( t/T ) * cos( pi * beta  / T * t ) /
//				( 1 - pow( 2*beta/T * t, 2 ) );
//	}
//	return { v/T, t1 };
//}

double root_raisedcosine( double t, double Fs, double beta )
{
	double T = 1 / Fs;
	if( t == 0 )
		return 1 - beta + 4 * beta / pi;
	else if( abs(t) == T/4/beta )
		return beta/sqrt(2) * ( (1+2/pi) * sin(pi/4/beta) +
								(1-2/pi) * cos(pi/4/beta) );
	else
		return ( sin( pi*t/T*(1-beta) )
			   + 4*beta*t/T * cos( pi*t/T*(1+beta) ) )
			   / ( pi*t/T * (1 - pow(4*beta*t/T,2)) );
}

rseq root_raisedcosine( int N, double ws, double beta )
{
	assert( abs(ws) <= pi );
	int n1    = -floor(N/2.);
	double Fs = ws/(2*pi);
	rvec v(N);
	for( int n = 0; n < N; n++ ) {
		double t = n + n1;
		v(n)     = root_raisedcosine( t, Fs, beta );
	}
	return { Fs*v, n1 };
}

//rseq root_raisedcosine( int N, double ws, double beta )
//{
//	int t1   = -floor(N/2.);
//	double T = 2*pi / ws;
//	rvec v(N);
//	for( int i = 0; i < N; i++ ) {
//		double t = i + t1;
//		if( t == 0 )
//			v(i) = 1 - beta + 4 * beta / pi;
//		else if( abs(t) == T/4/beta )
//			v(i) = beta/sqrt(2) * ( (1+2/pi) * sin(pi/4/beta) +
//					                (1-2/pi) * cos(pi/4/beta) );
//		else
//			v(i) = ( sin( pi*t/T*(1-beta) )
//			       + 4*beta*t/T * cos( pi*t/T*(1+beta) ) )
//			       / ( pi*t/T * (1 - pow(4*beta*t/T,2)) );
//	}
//	return { v/T, t1 };
//}

#ifndef BEALAB_NOPYTHON
rseq hamming( int N )
{
	python::function<rvec> hamming( "scipy", "hamming" );
	rvec x = hamming( N );
	int t1 = -(N-1)/2;
	return rseq( x, t1 );
}
#endif

}
}
