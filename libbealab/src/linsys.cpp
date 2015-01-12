#include <bealab/extensions/control/linsys.hpp>

namespace bealab
{
namespace control
{

state_space discretize_sampling( const state_space& ss, double T )
{
	int N   = ss.x.size();
	rmat Ad = eye(N) + T * ss.A;
	rmat Bd = T * ss.B;
	return state_space( Ad, Bd, ss.C, ss.D );
}

state_space controllable_canonical_form( const transfer_function& tf )
{
	// Check for properness
	rvec b = tf.num();
	rvec a = tf.den();
	if( b.size() > a.size() )
		error("controllable_canonical_form() - the transfer function is not proper");

	// Parse numerator and denominator
	rvec ar = a(range(1,a.size()));
	int N   = ar.size();
	b       = { zeros(N+1-b.size()), b };
	rvec br = b(range(1,b.size()));

	// Matrix A
	rmat A = zeros(N,N);
	A.row(0) = -ar;
	A( range(1,N), range(0,N-1) ) = eye(N-1);

	// Matrix B
	rmat B = zeros(N,1);
	B(0,0) = 1;

	// Matrix C
	rmat C(1,N);
	C.row(0) = br;

	// Matrix D
	rmat D(1,1);
	D(0,0) = b(0);

	return {A,B,C,D};
}

//------------------------------------------------------------------------------
// Array support functions
//------------------------------------------------------------------------------
Mat<transfer_function> tfm2mtf( const arma<rmat,rvec> &S )
{
	// Take numerator and denominator
	Vec<rmat> b = S.num();
	Vec<rmat> a = S.den();

	// Constants
	int I  = b(0).size1();
	int J  = b(0).size2();

	// Transform numerator
	Mat<rvec> mb(I,J);
	int Kb = b.size();
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ ) {
			mb(i,j).resize(Kb);
			for( int k = 0; k < Kb; k++ )
				mb(i,j)(k) = b(k)(i,j);
		}

	// Transform denominator
	Mat<rvec> ma(I,J);
	int Ka = a.size();
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < I; j++ ) {
			ma(i,j).resize(Ka);
			for( int k = 0; k < Ka; k++ )
				ma(i,j)(k) = a(k)(i,j);
		}

	// Form the matrix of ARMA models and return it
	Mat<transfer_function> M(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ ) {
			int l = ma(i,j).size();
			M(i,j).set_coeffs( mb(i,j), {{1},ma(i,j)(range(1,l))} );
		}
	return M;
}

arma<rmat,rvec> mtf2tfm( const Mat<transfer_function> &M )
{
	// Constants
	int I = M.size1();
	int J = M.size2();
	imat KKb(I,J);
	imat KKa(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ ) {
			KKb(i,j) = M(i,j).num().size();
			KKa(i,j) = M(i,j).den().size();
		}
	int Kb = max(KKb);
	int Ka = max(KKa);

	// Transform numerator
	Vec<rmat> B(Kb);
	for( int k = 0; k < Kb; k++ ) {
		B(k).resize(I,J);
		for( int i = 0; i < I; i++ )
			for( int j = 0; j < J; j++ )
				B(k)(i,j) = M(i,j).num()(k);
	}

	// Transform denominator
	Vec<rmat> A(Ka);
	for( int k = 0; k < Ka; k++ ) {
		A(k).resize(I,I);
		for( int i = 0; i < I; i++ )
			for( int j = 0; j < I; j++ )
				A(k)(i,j) = M(i,j).den()(k);
	}

	// Form the matrix ARMA model and return it
	arma<rmat,rvec> S;
	if( Ka > 1 )
		S.set_coeffs( B, { {eye(I)}, A(range(1,A.size())) } );
	else
		S.set_coeffs( B, {eye(I)} );
	return S;
}

}
}
