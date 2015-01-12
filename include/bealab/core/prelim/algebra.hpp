/// @file bealab/core/prelim/algebra.hpp
/// Operation for all basic types used to build algebras.

#ifndef _BEALAB_PRELIM_ALGEBRA_
#define	_BEALAB_PRELIM_ALGEBRA_

namespace bealab
{
/// @defgroup prelim_algebra Abstract algebra
/// Operations to build algebraic structures using double and complex types.
/// @{

/// @name Real numbers (double)

/// Adjoint operation
inline
double adjoint( double x ) { return x; }

/// Inverse
inline
double inv( double x ) { return 1/x; }

/// Norm
inline
double norm( double x, int p=2 ) { return p!=0 ? abs(x) : (x==0 ? 0 : 1); }

/// Inner product
inline
double inner_prod( double x, double y ) { return x*y; }

/// Outer product
inline
double outer_prod( double x, double y ) { return x*y; }
/// @}

/// @name Complex numbers

/// Adjoint operation
inline
complex adjoint( const complex &x ) { return conj(x); }

/// Inverse
inline
complex inv( const complex &x ) { return 1/x; }

/// Norm
inline
double norm( const complex &x, int p=2 ) { return p!=0 ? abs(x) : (abs(x)==0 ? 0 : 1); }

/// Inner product
inline
complex inner_prod( const complex &x, const complex &y ) { return x*conj(y); }

/// Outer product
inline
complex outer_prod( const complex &x, const complex &y ) { return x*conj(y); }
/// @}

/// @}
}
#endif
