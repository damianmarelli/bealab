// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

#include <bealab/scilib/stats.hpp>
#include <bealab/scilib/fourier.hpp>
#include <bealab/extensions/signal/tfanalysis.hpp>
#include <bealab/extensions/signal/misc.hpp>

namespace bealab
{
namespace signal
{

#ifndef BEALAB_NOPYTHON
cmat spectrogram( rseq x, int M, int D, int N )
{
	Vec<rseq> X = fb2pp_signal( x, D );
	rseq h0     = hamming( N );
	Vec<cseq> h = filterbank_dtft( h0, M );
	Mat<cseq> H = fb2pp_filterbank( h, D );
	Vec<cseq> Y = H * X;
	int K  = max( entrywise( [](const cseq& x) { return x.size(); } )(Y) );
	int t1 = min( entrywise( [](const cseq& x) { return x.t1(); } )(Y) );
	cmat A(M,K);
	for( int m = 0; m < M; m++ )
		for( int k = 0; k < K; k++ )
			A(m,k) = Y(m)(k+t1);
	return A;
}
#endif

Vec<cseq> filterbank_dtft( const cseq& h0, int M )
{
	Vec<cseq> h(M);
	for( int m = 0; m < M; m++ )
		h(m) = frequency_shift( h0, 2*pi*m/(double)M );
	return h;
}

cseq quasidual_window( const cseq& h, int M, int D, int t1, int t2 )
{
	int mm = -floor( (h.t2()-t1) / (double)M );
	int mM = floor( (t2-h.t1()) / (double)M );

	cmat R = zeros( D*(mM-mm+1), t2-t1+1 );
	cvec d = zeros( D*(mM-mm+1) );

	int ir = 0;
	for( int im = 0; im < mM-mm+1; im++ ) {
		int m  = im + mm;
		for( int t = 0; t < D; t++ ) {
			const cseq& h1 = time_shift( h, m*M+t );
			const cseq& h2 = downsample( h1, D );
			const cseq& h3 = upsample( h2, D );
			cseq h4        = time_shift( h3, -t );
			h4             = h4.trunc( t1, t2 );
			cvec rw        = conj(h4).vec();
			if( norm(rw) > 0 ) {
				R.row(ir) = rw;
				d(ir) = (m==0)/(double)M;
				ir = ir+1;
			}
		}
	}

	R = R( range(0,ir), range(0,R.size2()) );
	d = d( range(0,ir) );

	const cvec& f_ = lls( R, d );
	cout << "noncannonical dual: err = " << norm(d-R*f_) << endl;
	return cseq( f_, t1 );
}

tuple<double,double> frame_bounds( const cseq& h0, int M, int D )
{
	// Constants
	int N = h0.size();

	// Compute U(d) = \sum_m |h_m(z)| |h_m(\Omega^d z)|
	Vec<rvec> U(D);
	rvec Beta(D);
	for( int d = 0; d < D; d++ ) {
		const cseq& h00 = h0;
		const cseq& H00 = abs( dtft( h00, N ) );
		const cseq& H0d = circshift( H00.vec(), round((double)N/D)*d );
		const rvec& U0  = real( element_prod( H00.vec(), H0d.vec() ) );
		U(d) = U0;
		for( int m = 1; m < M; m++ )
			U(d) += circshift( U0, round((double)N/M)*m );
		Beta(d) = max( U(d) );
	}

	// Compute the intermediate terms
	double A_ = min( U(0) );
	double B_ = max( U(0) );
	double R_ = 0;
	for( int d = 1; d < D; d++ )
		R_ += sqrt( Beta(d) * Beta(-d) );

	// Compute the bounds
	double A = (A_ - R_) / D;
	double B = (B_ + R_) / D;

	return tuple<double,double>( A, B );
}

double dual_window_error( const cseq& h0, const cseq& f0, int M, int D )
{
	const rseq& x      = randn(10e3);
	const Vec<cseq>& X = analysis( x, h0, M, D );
	const rseq& y      = real( synthesis( X, f0, D ) );
	return norm( x - y ) / norm( x );
}

cseq connonical_dual_window( const cseq& h0, int M, int D, double tol )
{
	auto bounds        = frame_bounds( h0, M, D );
	double A           = get<0>(bounds);
	double B           = get<1>(bounds);
	double K           = 2. / (A+B);
	cseq s             = h0;
	cseq f0;
	for(;;) {
		f0 += K * s;
		double err = dual_window_error( h0, f0, M, D );
//		cout << "cannonical_dual() - Error = " << err << endl;
		if( err < tol )
			break;
		const Vec<cseq>& H = analysis( s, h0, M, D );
		s   = s - K * synthesis( H, h0, D );
	}
	f0 = f0.trim( .1*tol*max(abs(f0).vec()) );

	return f0;
}

}
}
