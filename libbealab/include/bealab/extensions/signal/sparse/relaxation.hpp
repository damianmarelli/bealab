/// @file bealab/extensions/signal/sparse/relaxation.hpp
/// Sparsify by solving a constrained optimization problem.

#ifndef _BEALAB_EXTENSIONS_SIGNAL_SPARSE_RELAXATION_
#define	_BEALAB_EXTENSIONS_SIGNAL_SPARSE_RELAXATION_

#include <bealab/extensions/signal/sparse/base.hpp>

namespace bealab
{
namespace sparse
{
/// @defgroup sparse-relaxation Relaxation algorithms
/// Sparsify by solving a constrained optimization problem.
/// @{

/// Base for all relaxation algorithms (abstract class).
class relaxation : public virtual base<double> {

	// Hide members
	using base::max_sparsity;

protected:

	/// Optimize the relaxed problem
	virtual
	rvec optimize() = 0;

public:

	double zero_threshold = 1e-6;												///< Threshold for zero-rounding the resulting coefficients
	ivec zth_indexes;															///< If non-empty, the zero-thresholding is done only on these indexes
	rvec guess;																	///< Initial parameter guess

	/// Destructor
	virtual
	~relaxation() {}

	/// Approximation error function
	function<double( const rvec& )> errfun;

	/// Gradient of the approximation error function
	function<rvec( const rvec& )> errfun_gradient = [this]( const rvec& coeffs )
	{
		return gradient( [this,&coeffs](const rvec& x){ return this->errfun(x); }, coeffs );
	};

	/// Virtual function override
	void approximate() override;
};

/// Generic lp relaxation sparse problem.
class lp_relaxation : public virtual relaxation {

	/// Direct approximatio algorithm.
	/// Minimize || x ||_p subject to errfun(x) < tol
	rvec optimize_direct();

	/// Recasted approximatio algorithm.
	/// Minimize || b ||_p subject to errfun(x) < tol and -b < x < b
	rvec optimize_recast();

	/// Virtual function override
	rvec optimize() override;

public:

	bool direct_method = false;													/// Flag to choose between the direct/recast method
	double p = 1;																///< p-norm to minimize
};

/// Vector/matrix version of lp-relaxation
class lp_relaxation_vector : public lp_relaxation {

	rmat gramian;																///< Gramian matrix
	rvec innerprods;															///< Inner products of the target with all atoms

public:

	rvec target;																///< Target
	rmat matrix;																///< Synthesis matrix

	/// Constructor
	lp_relaxation_vector();

	/// Virtual function override
	void approximate() override;
};

/// Hyder-Mahata method
class shrinking_gaussian : public virtual relaxation {

	/// Optimization for a given width
	rvec optimize( double width, const rvec& previous );

protected:

	/// Virtual function override
	rvec optimize() override;
};

///// Constrained sparsification combining shrinking Gaussian with the barrier method
//class shrinking_gaussian : public relaxation {
//
//	/// Gradient of the Gaussian weighting function
//	function<rvec(const rvec&,double)> grad_g = [this]( const rvec& x, double gw )
//	{
//		int N  = zth_indexes.size();
//		rvec g = zeros(x.size());
//		for( int n = 0; n < N; n++ ) {
//			int idx   = zth_indexes(n);
//			double xn = x(idx);
//			g(idx)    = xn / pow(gw,2) * exp( -pow(xn/gw,2) / 2 );
//		}
//		return g;
//	};
//
//	/// Graient of the barrier function
//	function<rvec(const rvec&)> grad_b = [this]( const rvec& x )
//	{
//		return errfun_gradient(x) / (tolerance - errfun(x));
//	};
//
//	/// Optimization for a given width
//	rvec optimize( double t, double gw, const rvec& wstart );
//
//	/// Virtual function override
//	rvec optimize() override;
//};

/// @}
}
}
#endif
