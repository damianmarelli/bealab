/// @file bealab/extensions/control/sysid.hpp
/// System identification library.

#ifndef _BEALAB_SYSID_
#define	_BEALAB_SYSID_

#include <bealab/extensions/control/linsys.hpp>

namespace bealab
{
namespace control
{
/// System identification module
namespace sysid
{
//------------------------------------------------------------------------------
/// @defgroup sysid System identification
/// System identification functions.
/// @{

/// Base class for system identification methods (virtual).
class sysid_b {
protected:

	int n_num;																	///< Number of numerator coefficients
	int	n_den;																	///< Number of denominator coefficients

public:

	rvec num;																	///< Identified numerator
	rvec den;																	///< Identified denominator
	int trace = 0;																///< Tracing level

	sysid_b( int nn, int nd ) { n_num = nn; n_den = nd; }						///< Constructor. Sets the numerator and denominator orders
	virtual ~sysid_b() {}														///< Destructor
	virtual void identify() = 0;												///< Identification function (virtual)
};

/// Base class for time-domain methods (virtual).
class timedomain_b : public sysid_b {
protected:
	rvec in;																	///< Input signal
	rvec out;																	///< Output signal
public:
	using sysid_b::sysid_b;
	void input( const rvec &in_ ) { in = in_; }									///< Set the input signal
	void output( const rvec &out_ ) { out = out_; }								///< Set the output signal
};

/// Time-domain identification using Kalman's method.
class kalman_td : public timedomain_b {
public:
	using timedomain_b::timedomain_b;
	void identify() override;													///< Virtual function override
};

/// Time-domain identification using Steiglitz's method.
class steiglitz_mcbride_td : public timedomain_b {
public:
	using timedomain_b::timedomain_b;
	void identify() override;													///< Virtual function override
};

/// Base class for frequency-domain methods (virtual).
class freqdomain_b : public sysid_b {
protected:
	cvec in;																	///< Input spectrum
	cvec out;																	///< Output spectrum
	rvec w;																		///< Input/output frequencies
public:
	using sysid_b::sysid_b;
	void input( const cvec &in_ ) { in = in_; }									///< Set the input spectrum
	void output( const cvec &out_ ) { out = out_; }								///< Set the output spectrum
	void frequencies( const rvec &w_ ) { w = w_; }								///< Set the input/output frequencies
};

/// Frequency-domain identification using Kalman's method.
class kalman_fd : public freqdomain_b {
public:
	using freqdomain_b::freqdomain_b;
	void identify() override;													///< Virtual function override
};

/// Frequency-domain identification using iterative reweighting.
class iterative_reweighting : public freqdomain_b {
	cvec apply_fun( const cvec&, const cvec& );									///< Apply fun() entry-wise
public:
	using freqdomain_b::freqdomain_b;
	function<complex(complex,complex)> fun;										///< Function whose square sum is to be minimized
	void identify() override;													///< Virtual function override
};

/// Frequency-domain identification using Steiglitz's method.
class steiglitz_mcbride_fd : public iterative_reweighting {
	using iterative_reweighting::fun;
public:
	steiglitz_mcbride_fd( int nn, int dd );										///< Constructor
};

/// Frequency-domain identification using Kobayashi's method.
class kobayashi : public iterative_reweighting {
	using iterative_reweighting::fun;
public:
	kobayashi( int nn, int dd );												///< Constructor
};

/// Frequency-domain identification using quasi-Newton.
class quasinewton : public freqdomain_b {
	double objfun(const rvec&);													///< Objective function to be minimized
	rvec objgrad(const rvec&);													///< Gradient of the objective function
public:
	using freqdomain_b::freqdomain_b;
	function<complex(complex,complex)> fun;										///< Function whose square sum is to be minimized
	function<complex(complex,complex)> dfun;									///< Derivative of fun on the second entry
	void identify() override;													///< Virtual function override
};

/// Frequency domain identification using quasi-Newton optimization and
/// linear least squares
class linear_fd : public quasinewton {
	using quasinewton::fun;
	using quasinewton::dfun;
public:
	linear_fd( int nn, int dd );												///< Constructor
};

/// Frequency domain identification using quasi-Newton optimization and
/// logarithmic least squares
class logarithmic_fd : public quasinewton {
	using quasinewton::fun;
	using quasinewton::dfun;
public:
	logarithmic_fd( int nn, int dd );											///< Constructor
};

/*
/// ARMA model identification using pseudo-linear regression
class plinreg : public arma_b {

	/// Internal iteration
	void iteration();

public:

	/// Virtual function override
	void identify() override;
};
*/

/// @}
}
}
}
#endif
