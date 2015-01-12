#include <fftw3.h>
#include <bealab/scilib/fourier.hpp>

namespace bealab
{

// Static mutex for calling fftw_plan_dft_1d() and fftw_destroy_plan()
static mutex dft_mutex;

Vec<complex> dft( const Vec<complex>& x_ )
{
	// Output
	cvec& x = const_cast<cvec&>( x_ );
	int N   = x.size();
	cvec y(N);

	// Input & output buffers
	fftw_complex* in  = reinterpret_cast<fftw_complex*>(&x(0));
	fftw_complex* out = reinterpret_cast<fftw_complex*>(&y(0));

	// Allocate a plan
	fftw_plan p;
	{
		lock_guard<mutex> lock(dft_mutex);
		p = fftw_plan_dft_1d( N, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
	}

	// FFT
	fftw_execute_dft( p, in, out );

	// Free the plan
	{
		lock_guard<mutex> lock(dft_mutex);
		fftw_destroy_plan(p);
	}

	return y;
}

Vec<complex> dft( const Vec<double>& x )
{
	cvec x_ = x;
	return dft( x_ );
}

Vec<complex> idft( const Vec<complex>& x_ )
{
	// Output
	cvec& x = const_cast<cvec&>( x_ );
	int N   = x.size();
	cvec y(N);

	// Input & output buffers
	fftw_complex *in  = reinterpret_cast<fftw_complex*>(&x(0));
	fftw_complex *out = reinterpret_cast<fftw_complex*>(&y(0));

	// Allocate a plan
	fftw_plan p;
	{
		lock_guard<mutex> lock(dft_mutex);
		p = fftw_plan_dft_1d( N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE );
	}

	// IFFT
	fftw_execute_dft( p, in, out );

	// Free the plan
	{
		lock_guard<mutex> lock(dft_mutex);
		fftw_destroy_plan(p);
	}

	return y / N;
}

Vec<complex> idft( const Vec<double>& x )
{
	cvec x_ = x;
	return idft( x_ );
}

rvec unwrap( rvec p )
{
	int N = p.size();
	for( int n = 0; n < N-1; n++ ) {
		if( p(n+1) - p(n) > pi )
			p( range( n+1, p.size() ) ) -= 2*pi * ones(p.size()-n-1);
		else if( p(n+1) - p(n) < -pi )
			p( range( n+1, p.size() ) ) += 2*pi * ones(p.size()-n-1);
	}
	return p;
}

cmat dft( int N )
{
	cmat W(N,N);
	for( int n = 0; n < N; n++ )
		W.column(n) = dft(rvec(unit(N,n)));
	return W;
}

cmat idft( int N )
{
	cmat W(N,N);
	for( int n = 0; n < N; n++ )
		W.column(n) = idft(rvec(unit(N,n)));
	return W;
}

}
