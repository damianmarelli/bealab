// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/extensions/signal/sparse/greedy.hpp
/// Choose the coefficients sequentially.

#ifndef _BEALAB_EXTENSIONS_SIGNAL_SPARSE_GREEDY_
#define	_BEALAB_EXTENSIONS_SIGNAL_SPARSE_GREEDY_

#include <bealab/extensions/signal/sparse/base.hpp>

namespace bealab
{
namespace signal
{
namespace sparse
{
/// @defgroup sparse-greedy Greedy algorithms
/// Choose the coefficients sequentially.
/// @{

/// Base for all greedy sparse algorithms (abstract class).
/// @param RAN range type
/// @param COE coefficient type
/// @param IDX index type
template<class COE=double>
class greedy : public virtual base<COE> {

	figure fig = figure("Error");												///< Figure for plotting

protected:

	/// Error function
	virtual
	double errfun( const ivec&, const vec<COE>& ) = 0;

	/// Function to force a stop of the iterations
	virtual
	bool forcestop() { return false; }

	/// Returns the next set of indexes
	virtual
	ivec next_indexes() = 0;

	/// Optimal coefficients for the given indexes
	virtual
	vec<COE> optimal_coefficients() = 0;

	/// Virtual function override
	void plotfun() override
	{
		fig.clear().logy().rangex(1,nan).plot( sparsity_history, error_history );
	}

public:

	///@name History
	vec<vec<COE>> coefficient_history;											///< Coefficient history
	rvec error_history;															///< Error history
	rvec sparsity_history;														///< Sparsity history
	///@}

	/// Virtual destructor
	virtual
	~greedy() {}

	/// Virtual function override
	void approximate() override
	{
		// Initialize the number of indexes
		this->indexes.resize(0);

		// Main loop
		error_history.resize(0);
		sparsity_history.resize(0);
		for( int i = 0;; i++ ) {

			// Get optimal coefficients
			this->coefficients = optimal_coefficients();
			if( i > 0 )
				coefficient_history = { coefficient_history,
						vec<rvec>{this->coefficients} };

			// Evaluate the approximation error
			this->error = errfun( this->indexes, this->coefficients );

			// Trace
			if(this->trace)
				cout << "Iteration: " << i << ", error: " << this->error << endl;

			// Plot
			error_history    = { error_history, {this->error} };
			sparsity_history = { sparsity_history, {(double)this->indexes.size()} };
			if( this->plot )
				this->plotfun();

			// Check stop condition
			if( this->error < this->tolerance ||
					(int)this->indexes.size() >= this->max_sparsity ||
					this->forcestop() )
				break;

			// Add one extra index
			this->indexes = ivec{ this->indexes, next_indexes() };
		}
	}
};

/// Generic orthogonal matching pursuit.
template<class COE=double, class RAN=vec<COE>>
class orthogonal_matching_pursuit : public virtual greedy<COE> {

	mat<COE> grammian;															///< Grammian matrix with the current atoms
	vec<COE> innerps;															///< Inner-products between the target and the current atoms
	vec<COE> all_innerps;														///< All inner-products between the target and the atoms
	rvec all_norms;																///< All atom norm

protected:

    /// Norm of a vector in the range space
	virtual
    double norm( const RAN& X )
    {
    	return bealab::norm( X );
    }

	/// Virtual function override
	double errfun( const ivec& idxs, const vec<COE>& coeffs ) override
	{
		approximation = synthesis( idxs, coeffs );
		return norm( approximation - target );
	}

	figure figc;

	/// Virtual function override
    ivec next_indexes() override
	{
		// Compute all analysis coefficients
		vec<COE> ainnerps = analysis( this->target - this->approximation );

		// Build the normalized correlations
		vec<COE> ncorrs = abs( ainnerps );
		int I = ncorrs.size();
		for( int i = 0; i < I; i++ )
			ncorrs(i) /= ( all_norms(i) < 1e-6 ? inf : all_norms(i) );

		// Pick best integer index
		int idx;
		for(;;) {

			// Maximum
			idx = max_index( ncorrs );

			// Make sure we don't choose an already chosen index
			auto it = find( this->indexes.begin(), this->indexes.end(), idx );
			if( it == this->indexes.end() )
				break;
			else {
				ncorrs(idx) = 0;
				if(this->trace)
					cout << "*** Rejected index: " << idx << endl;
			}
		}

		// Return the index
		return {idx};
	}

	/// Virtual function override
	vec<COE> optimal_coefficients() override
	{
		// Current number of indexes
		int I1 = this->coefficients.size();
		int I2 = this->indexes.size();
		grammian.resize(I2,I2);
		innerps.resize(I2);

		// Update Gramian matrix
		for( int i = I1; i < I2; i++ ) {

			for( int j = 0; j < i; j++ ) {
				grammian(j,i) = inner_product_atoms( this->indexes(i), this->indexes(j) );
				grammian(i,j) = inner_product_atoms( this->indexes(j), this->indexes(i) );
			}
			grammian(i,i) = inner_product_atoms( this->indexes(i), this->indexes(i) );

			// Update current inner-products between the target and atoms
			innerps(i) = all_innerps( this->indexes(i) );
		}

		// Return optimal sparse coefficients
		return linsolve( grammian, innerps );
//		return lls( grammian, innerps );
	}

public:

	/// Virtual destructor
	virtual
	~orthogonal_matching_pursuit() {}

    /// @name Range space
    RAN target;																	///< Target
	RAN approximation;															///< Final approximation
	/// @}

    /// @name Functors

	/// Computes the inner products between an element and all the atoms
    function<vec<COE>( const RAN& )> analysis;

	/// Build a sparse approximation
	function<RAN( const ivec&, const vec<COE>& )> synthesis;

	/// Inner product of atoms given the indexes
    function<COE( int, int )> inner_product_atoms;
	/// @}

	/// Initialization
	void initialization()
	{
		// Reset
		grammian.resize(0,0);
		innerps.resize(0);
		this->coefficients.resize(0);

		// Precompute all inner-products with the target
		all_innerps = analysis( this->target );

		// Precomputes the norm of all atoms
		int I = all_innerps.size();
		all_norms.resize(I);
		for( int i = 0; i < I; i++ )
			all_norms(i) = sqrt( abs(inner_product_atoms(i,i)) );
	}

    /// Virtual function override
	void approximate() override
	{
		initialization();
		greedy<COE>::approximate();
	}
};

/// Vector/matrix version of orthogonal matching pursuit
template<class COE=double>
class orthogonal_matching_pursuit_vector : public orthogonal_matching_pursuit<COE> {

	mat<COE> gramian;															///< Gramian matrix

public:

    rmat matrix;																///< Synthesis matrix

    /// Constructor
    orthogonal_matching_pursuit_vector()
    {
		this->analysis = [this]( const vec<COE>& y )
		{
			return adjoint(this->matrix) * y;
		};

		this->synthesis = [this]( const ivec& idxs, const vec<COE>& x )
		{
			return rmat(this->matrix( indirect(vrange(0,this->matrix.size1())), indirect(idxs) )) * x;
		};

		this->inner_product_atoms = [this]( int i, int j )
		{
			return this->gramian(j,i);
		};
    }

	/// Virtual function override
	void approximate() override
	{
		gramian = adjoint(matrix) * matrix;
		orthogonal_matching_pursuit<COE>::approximate();
	}
};

/// Generic non-linear pursuit.
template<class RAN=rvec>
class gradient_pursuit : public virtual greedy<double> {

	int N;																		///< Total number of coefficients

protected:

	enum gradient_type { full, reduced };

	/// Gradient of the cost function.
	/// The indexes of the given coefficients 'coeff' are this->indexes.
	/// If whole = true, it returns the gradient entries corresponding this->all_indexes.
	/// If whole = false, it returns those corresponding to this->indexes.
	virtual
	rvec gradient( const rvec& coeffs, gradient_type gt )
	{
		if( gt == reduced ) {

			// Objective function of the coefficients corresponding to either all_indexes or indexes
			function<double(const rvec&)> fun = [this]( const rvec& y ) -> double
			{
				return this->errfun( this->indexes, y );
			};

			// Return the gradient
			return bealab::gradient( fun, coeffs );
		}
		else {

			// Objective function of the coefficients corresponding to this->all_indexes
			function<double(const rvec&)> fun = [this]( const rvec& y ) -> double
			{
				return this->errfun( vrange(0,this->N), y );
			};

			// Build the vector at which the gradient will be evaluated
			rvec x = zeros( N );
			x(indirect(this->indexes)) = coeffs;

			// Return the gradient
			return bealab::gradient( fun, x );
		}
	}

	/// Virtual function override
	double errfun( const ivec& idxs, const rvec& coeffs ) override
	{
		return error_function( idxs, coeffs );
	}

	/// Virtual function override
	ivec next_indexes() override
	{
		// Compute the gradient
		rvec grad = gradient( this->coefficients, full );

		// Find maximum
		int idx;
		for(;;) {

			// Maximum
			idx = max_index( abs(grad) );

			// Make sure we don't choose an already chosen index
			auto it = find( this->indexes.begin(), this->indexes.end(), idx );
			if( it == this->indexes.end() )
				break;
			else {
				grad(idx) = 0;
				cout << "*** Rejected index: " << idx << endl;
			}
		}

		// Return the IDX index
		return ivec({idx});
	}

	/// Guess for the non-linear optimization
	virtual
	rvec guess()
	{
		vec<double> x0 = this->coefficients;
		int I = this->indexes.size() - this->coefficients.size();
		for( int i = 0; i < I; i++ )
			x0 = rvec{ x0, {0} };
		return x0;
	}

	/// Virtual function override
	rvec optimal_coefficients() override
	{
		// Deal with the case when there are zero indexes
		if( this->indexes.size() == 0 )
			return rvec();

		// Guess
		rvec x0 = guess();

		// Objective function
		function<double(const rvec&)> fun = [this]( const rvec& x ) -> double
		{
			return this->errfun( this->indexes, x );
		};

		// Derivative
		function< rvec(const rvec&)> der = [this]( const rvec& x ) -> rvec
		{
			return this->gradient( x, reduced );
		};

		// If the guess is not good, return it without optimizing.
		double testval = fun(x0);
		if( testval == inf || isnan(testval) )
			return x0;

		// Compute optimal coefficients
		optimization::quasinewton opt( x0.size() );
		opt.set_objective( fun, der );
		opt.guess   = x0;
		opt.trace   = trace - 1;
		rvec coeffs = opt.optimize();
//		if( opt.stopreason == opt.error ) {
//			cerr << "NLOPT failed!" << endl;
//			optimization::bfgs opt( x0.size() );
//			opt.set_objective( fun, der );
//			opt.guess = x0;
//			coeffs    = opt.optimize();
//		}
		return coeffs;
	}

public:

	/// Functor version of the error function
	function<double( const ivec&, const rvec& )> error_function;

	/// Set the number of coefficients
	void set_size( int N_ ) { N = N_; }
};

/// @}
}
}
}
#endif
