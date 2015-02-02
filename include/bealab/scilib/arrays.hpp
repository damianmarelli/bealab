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

/// @name Conversions mat/vec
template<class T>
mat<vec<T>> vm2mv( const vec<mat<T>>& x )
{
	if( x.size() == 0 )
		return mat<vec<T>>();
	mat<T> x0 = x(0);
	int I     = x0.size1();
	int J     = x0.size2();
	int N     = x.size();
	mat<vec<T>> y(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ ) {
			y(i,j).resize(N);
			for( int n = 0; n < N; n++ )
				y(i,j)(n) = x(n)(i,j);
		}
	return y;
}

template<class T>
vec<mat<T>> mv2vm( const mat<vec<T>>& x )
{
	int I = x.size1();
	int J = x.size2();
	if( x.size1() * x.size2() == 0 )
		return mat<vec<T>>(I,J);
	int N = x(0,0).zeros(I,J);
	vec<mat<T>> y(N);
	for( int n = 0; n < N; n++ ) {
		y(n).resize(I,J);
		for( int i = 0; i < I; i++ )
			for( int j = 0; j < J; j++ )
					y(n)(i,j) = x(i,j)(n);
	}
	return y;
}
/// @}

/// @name Conversions vec/sequence
template<class T>
vec<sequence<T>> sv2vs( const sequence<vec<T>> &x )
{
	vec<T> x1 = x(x.t1());
	int N = x1.size();
	vec<sequence<T>> y(N);
	for( int i = 0; i < N; i++ ) {
		auto component = [i](const vec<T>& x){return x(i);};
		y(i) = entrywise(component)( x );
	}
	return y;
}

template<class T>
sequence<vec<T>> vs2sv( const vec<sequence<T>> &x )
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

	sequence<vec<T>> y( t2-t1+1, t1 );
	for( int t = t1; t <= t2; t++ ) {
		auto sample = [t](const sequence<T>& xi){return xi(t);};
		y(t) = entrywise(sample)(x);
	}
	return y;
}
/// @}

/// @name Conversions mat/sequence
template<class T>
mat<sequence<T>> sm2ms( const sequence<mat<T>> &x )
{
	mat<T> x1 = x(x.t1());
	int I = x1.size1();
	int J = x1.size2();
	mat<sequence<T>> y(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ ) {
			auto component = [i,j](const mat<T>& x){return x(i,j);};
			y(i,j) = entrywise(component)( x );
		}
	return y;
}

template<class T>
sequence<mat<T>> ms2sm( const mat<sequence<T>> &x )
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

	sequence<mat<T>> y( t2-t1+1, t1 );
	for( int t = t1; t <= t2; t++ ) {
		auto sample_t = [t](const sequence<T>& xij){return xij(t);};
		y(t) = entrywise(sample_t)( x );
	}
	return y;
}
/// @}

/// @name Sampling
template<class T>
vec<T> sample_get( const vec<sequence<T> > &v, int t )
{
	int	I = v.size();
	vec<T>	w(I);
	for( int i = 0; i < I; i++ )
		w(i) = v(i)(t);
	return w;
}

template<class T>
void sample_set( vec<sequence<T> > &v, int t, const vec<T> &w )
{
	int	I = v.size();
	for( int i = 0; i < I; i++ )
		v(i)(t) = w(i);
}

template<class T>
mat<T> sample_get( const mat<sequence<T>> &m, int t )
{
	int	I = m.size1();
	int	J = m.size2();
	mat<T>	M(I,J);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < J; j++ )
			M(i,j) = m(i,j)(t);
	return M;
}

template<class T>
void sample_set( mat<sequence<T> > &m, int t, const mat<T> &M )
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
mat<sequence<T>> inv( const mat<sequence<T>>& G, double tol=1e-12 )
{
	int N = max( G.apply( [](const cseq& x){ return x.size(); } ) );			// Maximum sequence size in G
	mat<sequence<T>> Gi;
	for(;; N *= 2 ) {

		// Compute the inverse
		mat<cvec> Gf   = G.apply( [N](const cseq& x){ return dtft(x,N); } );
		vec<cmat> Gf_  = ms2sm( mat<cseq>(Gf) );
		vec<cmat> Gfi_ = element_inv( Gf_ );
		mat<cvec> Gfi  = sm2ms( sequence<cmat>(Gfi_) );
		Gi             = apply( [N](const cvec& x){ return idtft(x); } )( Gfi );

		// Evaluate
		double err = norm( G - G*Gi*G ) / norm( G );
		cout << "inv(mat<cseq>) - error = " << err << endl;
		if( err < tol )
			break;
	}

	return Gi;
}
/// @}

/// @}
}
#endif
