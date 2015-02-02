/// @file bealab/extensions/control/baytrack.hpp
/// Implementation of Bayesian tracking techniques.

#ifndef _BEALAB_BAYTRACK_
#define	_BEALAB_BAYTRACK_

#include <bealab/core/blas.hpp>
#include <bealab/scilib/stats.hpp>
#include <bealab/scilib/optimization.hpp>
#include <bealab/extensions/control/linsys.hpp>
#include <boost/circular_buffer.hpp>

namespace bealab
{
namespace control
{
//------------------------------------------------------------------------------
/// @defgroup baytrack Bayesian tracking
/// Implementation of Bayesian tracking techniques.
/// @{

using boost::circular_buffer;

/// Kalman filter
class kalman_filter : public state_space {
public:

	rvec x;																		///< Initial state mean
	rmat P;																		///< Initial state covariance

	/// Constructor
	kalman_filter( const rmat& A, const rmat& B, const rmat& C, const rmat& D,
			const rmat& Q, const rmat& R, const rmat& P_ = rmat() ) :
				state_space(A,B,C,D,Q,R)
	{
		int N = A.size1();
		x     = zeros(N);
		if( P_.size1() == 0 )
			P = zeros(N,N);
		else
			P = P_;
	}

	/// Constructor
	kalman_filter( const state_space& ss, const rmat& P = rmat() ) :
		kalman_filter( ss.A, ss.B, ss.C, ss.D, ss.Q, ss.R, P) {}

	/// Prediction step
	rvec prediction( const rvec& u )
	{
		x = A * x + B * u;
		P = noproxy(A * P) * trans(A) + Q;
		return x;
	}

	/// Update step
	rvec update( const rvec& y )
	{
		int N  = A.size1();
		rmat K = trans( linsolve( noproxy(C * P) * trans(C) + R, C * P ) );
		x      = x + K * ( y - C * x );
		P      = (eye(N) - K * C) * P;
		return x;
	}

	/// Filter one pair of input and output samples
	rvec operator()( const rvec& u, const rvec& y )
	{
		prediction( u );
		update( y );
		return x;
	}

	/// Filter an output sample, assuming zero input.
	rvec operator()( const rvec& y )
	{
		int M  = B.size2();
		rvec u = zeros(M);
		return (*this)( u, y );
	}

	/// Filter one pair of sequences of input and output samples
	vec<rvec> operator()( const vec<rvec>& U, const vec<rvec>& Y )
	{
		assert( U.size() == Y.size() );
		int T = U.size();
		vec<rvec> Xh(T);
		for( int t = 0; t < T; t++ )
			Xh(t) = (*this)( U(t), Y(t) );
		return Xh;
	}

	/// Filter a sequences of output samples, assuming zero input.
	vec<rvec> operator()( const vec<rvec>& Y )
	{
		int T = Y.size();
		vec<rvec> Xh(T);
		for( int t = 0; t < T; t++ )
			Xh(t) = (*this)( Y(t) );
		return Xh;
	}
};

/// Base class for implementing a particle filter.
/// The virtual member functions  need to be implemented in a derived class.
/// Need to implement either state_transition() or particle_prediction(), and
/// either likelihood_function() or particle_update().
template< class state_t = rvec, class obs_t = state_t, class out_t = obs_t >
class particle_filter_b {
protected:

	/// A particles path with it's associated weight
	struct particle_path {
		circular_buffer<state_t> trajectory;
		double weight;
	};

	rvec old_weights;															///< Weights before the last update step
	bool init_f = false;														///< Flag to indicate that the filter was initialized
	int time    = 0;															///< The current time used when running the filter

	/// A particles with it's associated weight
	struct particle {
		state_t state;
		double weight;
	};

	/// Draw a sample from the initial distribution
	virtual
	state_t draw_initial_sample() = 0;

	/// State transition function
	virtual
	state_t state_transition( const state_t& state ) { return state; }

	/// Likelihood of the state given the observation
	virtual
	double likelihood_function( const obs_t& obs, const state_t& state ) { return 1; }

	/// Default particle prediction function. Uses state_transition().
	virtual
	particle particle_prediction( const particle& p )
	{
		particle pp = p;
		pp.state    = state_transition( p.state );
		return pp;
	}

	/// Default particle update function. Uses likelihood_function().
	virtual
	particle particle_update( const obs_t& obs, const particle& p )
	{
		particle pu = p;
		pu.weight  *= likelihood_function( obs, p.state );
		return pu;
	}

	/// Output function
	virtual
	out_t output_function( const state_t& state ) = 0;

public:

	vector<particle_path> x;													///< Particles of the current predicted or updated state
	const int I;																///< Number of particles
	const int lag;																///< Amount of lag used for smoothing.
	double resample_period = 1;

	/// Constructor
	particle_filter_b( int I_, int lag_=0 ) :
		old_weights(I_), x(I_), I(I_), lag(lag_)
	{
		assert( I > 0 );
		for( int i = 0; i < I; i++ )
			x[i].trajectory = circular_buffer<state_t>(lag+1);
	}

	/// Destructor
	virtual
	~particle_filter_b() {}

	/// Implements the prediction step
	virtual
	void prediction()
	{
		#pragma omp parallel for schedule(dynamic)
		for( int i = 0; i < I; i++ ) {
			particle p  = { x[i].trajectory.back(), x[i].weight };
			particle pp = particle_prediction( p );
			x[i].trajectory.push_back( pp.state );
			x[i].weight = pp.weight;
		}
	}

	/// Implements the update step
	virtual
	void update( const obs_t& z )
	{
		#pragma omp parallel for schedule(dynamic)
		for( int i = 0; i < I; i++ ) {
			particle p             = { x[i].trajectory.back(), x[i].weight };
			particle pu            = particle_update( z, p );
			x[i].trajectory.back() = pu.state;
			old_weights(i)         = x[i].weight;
			x[i].weight            = pu.weight;
		}
		normalize();
	}

	/// Implements the normalization step
	virtual
	void normalize()
	{
		double K = 0;
		for( int i = 0; i < I; i++ )
			K += x[i].weight;
//		if( K == 0 ) {
//			warning("Ran out of particles");
//			for( int i = 0; i < I; i++ )
//				x[i].weight = old_weights(i);
//		}
//		else {
//			for( int i = 0; i < I; i++ )
//				x[i].weight /= K;
//		}
		for( int i = 0; i < I; i++ )
			x[i].weight /= K;
		if( K == 0 )
			warning("Ran out of particles");
	}

	/// Implements the re-sampling step
	virtual
	void resample()
	{
		// Store current particles
		auto x0 = x;

		// Build the distribution
		rvec weights(I);
		for( int i = 0; i < I; i++ )
			weights(i) = x[i].weight;
		std::discrete_distribution<> d( weights.begin(), weights.end() );

		// Re-sample
		for( int i = 0; i < I; i++ ) {
			int idx = d(_stats::engine);
			x[i].trajectory = x0[idx].trajectory;
			x[i].weight     = 1./I;
		}
	}

	/// Generates the output
	virtual
	out_t output()
	{
		vec<out_t> o(I);
		#pragma omp parallel for schedule(dynamic)
		for( int i = 0; i < I; i++ )
			o(i) = output_function( x[i].trajectory.front() ) * x[i].weight;
		return sum(o);
	}

	/// Draw the particles of the initial state
	void initialize()
	{
		init_f = true;
		for( int i = 0; i < I; i++ ) {
			x[i].trajectory.push_back( draw_initial_sample() );
			x[i].weight = 1./I;
		}
	}

	/// Filter one sample
	out_t operator()( const obs_t& z )
	{
		if( !init_f )
			initialize();
		update( z );
		out_t o = output();
		if( mod( time+1, resample_period ) == 0 )
			resample();
		prediction();
		time++;
		return o;
	}

	/// Filter a vector of samples
	vec<out_t> operator()( const vec<obs_t>& Z )
	{
		int T = Z.size();
		vec<out_t> Yh(T);
		for( int t = 0; t < T; t++ )
			Yh(t) = (*this)( Z(t) );
		return Yh;
	}
};

/// Particle approximation to a Kalman filter (for testing purposes)
class particle_kalman_filter : public particle_filter_b<rvec,rvec,rvec>, public kalman_filter {

	int N;
	rmat Qh, Rh, Ph;

	rvec draw_initial_sample() override
	{
		return Ph * randn(N);
	}

	rvec state_transition( const rvec& x ) override
	{
		return A * x + Qh * randn(N);
	}

	double likelihood_function( const rvec& y, const rvec& x ) override
	{
		return pdf( multivariate_normal( C*x, Rh ), y );
	}

	rvec output_function( const rvec& x ) override
	{
		return C*x;
	}

public:

	using particle_filter_b::prediction;
	using particle_filter_b::update;
	using particle_filter_b::operator();

	particle_kalman_filter( int I, const rmat& A, const rmat& B, const rmat& C, const rmat& D,
			const rmat& Q, const rmat& R, const rmat& P_ = rmat() ) :
				particle_filter_b(I),
				kalman_filter(A,B,C,D,Q,R,P_)
	{
		N  = A.size1();
		Qh = real(msqrt(Q));
		Rh = real(msqrt(R));
		Ph = real(msqrt(P));
	}

	particle_kalman_filter( int I, const state_space& ss, const rmat& P = rmat() ) :
		particle_kalman_filter(I,ss.A,ss.B,ss.C,ss.D,ss.Q,ss.R,P)
	{}
};

/// State associated to a Rao-Blackwellized particle filter.
struct RB_state_t {
	rvec position;
	rvec mean;
	rmat covariance;
};

/// Base class for implementing Rao-Blackwellized particle filters.
/// The virtual member functions  need to be implemented in a derived class.
class RB_particle_filter_b : public particle_filter_b< RB_state_t, rvec, rvec > {

	/// Non-linear state transition function
	virtual
	rvec f( const rvec& x ) = 0;

	/// Linear state transition function
	virtual
	rvec g( const rvec& x ) = 0;

	/// Measurement non-linear component
	virtual
	rvec h( const rvec& x ) = 0;

	/// Output non-linear component
	virtual
	rvec o( const rvec& x ) = 0;

	/// Non-linear state transition matrix
	virtual
	rmat F( const rvec& x ) = 0;

	/// Linear state transition matrix
	virtual
	rmat G( const rvec& x ) = 0;

	/// Measurement linear component
	virtual
	rmat H( const rvec& x ) = 0;

	/// Output linear component
	virtual
	rmat O( const rvec& x ) = 0;

	/// Non-linear process noise covariance
	virtual
	rmat U( const rvec& x ) = 0;

	/// Linear process noise covariance
	virtual
	rmat V( const rvec& x ) = 0;

	/// Output noise covariance
	virtual
	rmat W( const rvec& x ) = 0;

	/// Override
	particle particle_prediction( const particle& p ) override
	{
		// Predicted state
		particle pp;

		// Non-linear elements
		rvec f_ = f( p.state.position );
		rmat F_ = F( p.state.position );
		rmat U_ = U( p.state.position );
		rvec g_ = g( p.state.position );
		rmat G_ = G( p.state.position );
		rmat V_ = V( p.state.position );

		// Non-linear prediction
		rvec nlp_mean   = f_ + F_ * p.state.mean;
		rmat nlp_cov    = noproxy(F_ * p.state.covariance) * trans(F_) + U_;
		rmat nlp_cov_h  = real(msqrt(nlp_cov));
		rvec nlp        = rand( multivariate_normal( nlp_mean, nlp_cov_h ) );
		pp.state.position = nlp;

		// Linear prediction
		if( det(nlp_cov) != 0 ) {

			// Linear false update
			rvec lfu_obs  = nlp - f_;
			rmat K        = trans( linsolve( nlp_cov, F_ * p.state.covariance ) );
			rvec lfu_mean = p.state.mean + K * ( lfu_obs - F_ * p.state.mean );
			rmat lfu_cov  = (eye(K.size1()) - K * F_) * p.state.covariance;

			// Prediction
			pp.state.mean       = g_ + G_ * lfu_mean;
			pp.state.covariance = noproxy(G_ * lfu_cov) * trans(G_) + V_;
		}
		else {
			pp.state.mean       = g_ + G_ * p.state.mean;
			pp.state.covariance = noproxy(G_ * p.state.covariance) * trans(G_) + V_;
		}

		// Weight
		pp.weight = p.weight;

		return pp;
	}

	/// Override
	particle particle_update( const rvec& obs, const particle& p ) override
	{
		// Update particles
		particle pu;

		// Non-linear elements
		rvec h_ = h( p.state.position );
		rmat H_ = H( p.state.position );
		rmat W_ = W( p.state.position );

		// Linear output prediction
		rvec lop_mean = H_ * p.state.mean + h_;
		rmat lop_cov  = noproxy(H_ * p.state.covariance) * trans(H_) + W_;

		// Check the singularity of lop_cov
		double singular = (det(lop_cov) == 0);

		// Non-linear state update
		pu.state.position = p.state.position;
		rvec delta        = obs - lop_mean;
		if( !singular ) {
			int N     = lop_mean.size();
			double a  = inner_prod( delta, linsolve( lop_cov, delta ) );
			double b  = log( pow(2*pi,N) * det(lop_cov) );
			pu.weight = exp( -0.5 * (a+b) ) * p.weight;
			if( isnan(pu.weight) )
				error("Weight is nan");
		}
		else
			pu.weight = ( norm(delta) == 0 ? 1 : 0) * p.weight;

		// Linear state update
		if( !singular ) {
			rvec lsu_obs        = obs - h_;
			rmat K              = trans( linsolve( lop_cov, H_ * p.state.covariance ) );
			pu.state.mean       = p.state.mean + K * ( lsu_obs - H_ * p.state.mean );
			pu.state.covariance = (eye(K.size1()) - K * H_) * p.state.covariance;
		}
		else {
			pu.state.mean       = p.state.mean;
			pu.state.covariance = p.state.covariance;
		}

		return pu;
	}

	/// Override
	rvec output_function( const RB_state_t& state ) override
	{
		rvec o_ = o( state.position );
		rmat O_ = O( state.position );
		return O_ * state.mean + o_;
	}

public:

	/// Constructor
	RB_particle_filter_b( int I, int lag=0 ) :
		particle_filter_b<RB_state_t,rvec,rvec>( I, lag ) {}

	/// Destructor
	virtual
	~RB_particle_filter_b() {}
};

/// Base class for implementing a particle smoother.
/// The virtual member functions  need to be implemented in a derived class.
//template< class pos_t = rvec >
//class particle_smoother_b : public particle_filter_b<pos_t> {
//
//	int L;																		///< Lag size
//	particles current;															///< Current updated particles
//	deque<particles> lagbuffer;													///< Buffer with the updated particles within the lag
//
//	/// Run one lag step
//	particles smooth_one_step( const particles& filtered, const particles& lagged )
//	{
////		particles rv = filtered;
////		for( int i = 0; i < I; i++ )
////			rv[i].weight *=
//		return lagged;
//	}
//
//	/// Run the backwards smoothing process
//	particles smooth_backwards()
//	{
//		particles lagged = current;
//		for( int l = 0; l < L; l++ )
//			lagged = smooth_one_step( lagbuffer[l], lagged );
//		return lagged;
//	}
//
//public:
//
//	/// Constructor
//	particle_smoother_b( int I_, int L_ ) :
//		particle_filter_b(I_), L(L_) {};
//
//	/// Filter one sample
//	pos_t operator()( const rvec& z ) override
//	{
//		// Initialize
//		if( !init_f )
//			initialize();
//
//		// Update buffer
//		if( (int)lagbuffer.size() == L )
//			lagbuffer.pop_back();
//		lagbuffer.emplace_front( current );
//
//		// Filter
//		current             = update( x, z );
//		particles resampled = resample( current );
//		x                   = prediction( resampled );
//		time++;
//
//		// Smoothing
//		particles xs = smooth_backwards();
//		return mean( output(xs) );
//	}
//};

/// ML Kalman filter
class maximum_likelihood_kalman_filter_b {

	/// State type
	struct state {
		rvec x;
		rmat P;
	};

	rmat C;
	rmat C_pseudoinverse;

	virtual
	double logLF( const rvec& z, const rvec& x ) = 0;

	virtual
	rvec logLF_gradient( const rvec& z, const rvec& x )
	{
		auto loglf = [this,&z](const rvec&x){return logLF(z,x);};
		return gradient( loglf, x );
	}

	virtual
	rvec output( const rvec& x )
	{
		return x;
	}

	/// Maximum likelihood estimation.
	/// Returns the ML estimate in rv.x and the Hessian of the logLF in rv.P
	state ml_estimate( const rvec& z, const rvec& guess=rvec() )
	{
		// Detect Rao-Blackwellization
		int L = s.x.size();
		int N = s.x.size();
		if( C_pseudoinverse.size1() * C_pseudoinverse.size2() > 0 )
			N = C_pseudoinverse.size2();

		// ML Estimation
		auto objfun  = [this,&z](const rvec&x){return -logLF(z,x);};
		auto objgrad = [this,&z](const rvec&x){return -logLF_gradient(z,x);};
		optimization::bfgs opt(N);
		opt.set_objective( objfun, objgrad );
		if( guess.size() == 0 )
			opt.guess = zeros(N);
		else
			opt.guess = guess;

//		// Test derivatives
//		rvec xt = randn(N);
//		cout << "xtest = " << xt << endl;
//		opt.test_derivatives( xt );
//		cin.get();

//		opt.stop_fincrement_relative = 1e-2;
//		opt.stop_xincrement_relative = 1e-2;
		rvec xh = opt.optimize();

		// Return value
		state ml;
		ml.x = { xh, zeros(L-N) };
		ml.P = jacobian( objgrad, xh );
		ml.P = ( ml.P + trans(ml.P) ) / 2;
		ml.P = { { ml.P,         zeros(N,L-N)   },
				 { zeros(L-N,N), zeros(L-N,L-N) } };
		return ml;
	}

	/// Prediction step
	state prediction( const state& s )
	{
		state p;
		p.x = A * s.x;
		p.P = noproxy(A * s.P) * trans(A) + Q;
		return p;
	}

	/// Update step
	state update( const state& s, const rvec& z )
	{
		// ML estimate
		state ml = ml_estimate( z, C*s.x );

		// Update the state
		state u;
		u.P = inv( inv(s.P) + ml.P );
		u.x = u.P * ( linsolve( s.P, s.x ) + ml.P * ml.x );
		return u;
	}

public:

	state s;																	///< Current state
	rmat A;																		///< State transition matrix
	rmat Q;																		///< Process noise covariance

	/// Destructor
	virtual
	~maximum_likelihood_kalman_filter_b() {}

	void set_C( const rmat& C_ )
	{
		C = C_;
		C_pseudoinverse = real(pinv(C));
	}

	/// Filter the sample z
	rvec operator()( const rvec& z )
	{
		s = prediction( s );
		s = update( s, z );
		return output(s.x);
	}

	/// Filter a vector of samples
	vec<rvec> operator()( const vec<rvec>& Z )
	{
		int T = Z.size();
		vec<rvec> Yh(T);
		for( int t = 0; t < T; t++ )
			Yh(t) = (*this)( Z(t) );
		return Yh;
	}
};

/// @}
}
}
#endif
