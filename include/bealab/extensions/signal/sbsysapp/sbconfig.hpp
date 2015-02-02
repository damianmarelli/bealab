/// @file bealab/extensions/signal/sbsysapp/sbconfig.hpp
/// Models a subband index, a subband model and a subband configuration.

#ifndef _BEALAB_EXTENSIONS_SIGNAL_SBSYSAPP_SBCONFIG_
#define	_BEALAB_EXTENSIONS_SIGNAL_SBSYSAPP_SBCONFIG_

#include <bealab/core/matfile.hpp>
#include <bealab/scilib/sequences.hpp>
#include <bealab/scilib/arrays.hpp>
#include <bealab/extensions/control/linsys.hpp>
#include <bealab/extensions/signal/tfanalysis.hpp>

namespace bealab
{
namespace signal
{
namespace sbsysapp
{
/// @defgroup sbsysapp_sbconfig Subband modeling
/// Models a subband index, a subband model and a subband configuration.
/// @{

/// Subband index (m,n,t,real_f,num_f).
struct sbindex {
private:

	friend class boost::serialization::access;

	/// Serialization
	template<class Archive>
	void serialize( Archive& ar, const unsigned int version )
	{
		ar & m;
		ar & n;
		ar & t;
		ar & real_f;
		ar & num_f;
	}

public:
	int m;																		///< Row
	int n;																		///< Column
	int t;																		///< Time
	bool real_f;																///< Flag indicating if the index is real
	bool num_f;																	///< Flag indicating if the index belongs to the numerator

	/// Check if two indexes are the same
	bool operator==( const sbindex& idx ) const;

	/// Check if two indexes are different
	bool operator!=( const sbindex& idx ) const;

	/// Index comparison for sorted containers
	bool operator<( const sbindex& idx ) const;

	/// Get the conjugate of an index w.r.t. M
	sbindex conjugate( int M ) const;

	/// Check if an index is self-conjugate w.r.t. M
	bool is_selfconjugate( int M ) const;
};

/// Display a subband index
template<class A, class B>
std::basic_ostream<A,B> &operator<<( std::basic_ostream<A,B> &os, const sbsysapp::sbindex &sbi )
{
	string realimag = (sbi.real_f ? "real" : "imag");
	string numden   = (sbi.num_f  ? "num"  : "den");
	os << "{ " << sbi.m << ", " << sbi.n << ", " << sbi.t << ", " << realimag << ", " << numden << " }";
	return os;
}

/// A set of subband indexes with an associated number of subbands
class sbindex_pool {
protected:

	friend class boost::serialization::access;

	/// Serialization
	template<class Archive>
	void serialize( Archive& ar, const unsigned int version )
	{
		ar & sbindexes;
		ar & M;
	}

	bool real_target_f;															///< Flag to indicate that the target is real

	/// Convert a subband index into an integer index
	int sbi2int( const sbindex& idx );

public:

	vec<sbindex> sbindexes;														///< Pool of all subband indexes
	int M;																		///< Number of subbands

	/// Constructor
	sbindex_pool( int M_, bool rt=true ) :
		real_target_f(rt), M(M_) {}

	/// Virtual destructor
	virtual
	~sbindex_pool() {}

	/// Initialize table of subband indexes with an array of indexes
	virtual
	void set_sbindexes( const vec<sbindex>& idxs );

	/// Initialize table of subband indexes with a pattern matrix
	void set_sbindexes( const mat<cseq>& MA, const mat<cseq>& AR );

	/// Initialize table of subband indexes specifying number of off-diagonal terms
	void set_sbindexes( int noffdiags, int t1, int t2, bool den_f=false );

	/// Sparsity of the index pool
	int sparsity() const;
};

/// Subband ARMA model
class sbmodel : public sbindex_pool {

	friend class boost::serialization::access;
	friend struct matvar<sbmodel>;

	/// Serialization
	template<class Archive>
	void serialize( Archive& ar, const unsigned int version )
	{
        ar & boost::serialization::base_object<sbindex_pool>(*this);
		ar & sbvalues;
	}

	/// Model matrices
	struct models_t {
		mat<cseq> S;															///< MA part
		mat<cseq> P;															///< AR part
	};

	/// Computes the impulse response length of inv(I-P)
	int inv_ImP_irlength( const vec<cmat>& P ) const;

	/// Computes the impulse response of inv(I-P)
	mat<cseq> inv_ImP( const mat<cseq>& P_ ) const;

public:

	rvec sbvalues;																///< Values associated to the subband indexes

	/// Constructor
	sbmodel( int M_=0 ) : sbindex_pool(M_) {}

	/// Obtain the subband model matrices
	models_t models() const;

	/// Obtain the impulse response
	mat<cseq> impulse_response() const;
};

/// A subband configuration (with a single subband model).
/// It is formed by a subband model, analysis and synthesis filterbanks
/// and the downsampling factor.
class sbconfig {

	cseq _h0;																	///< Analysis prototype
	cseq _f0;																	///< Synthesis prototype

protected:

	static const int callback_trace_level = 100;
	mat<cseq> H;																///< Polyphase matrix of the analysis filterbank
	mat<cseq> F;																///< Polyphase matrix of the synthesis filterbank
	mat<cseq> FA;																///< Pre-computed F.A() to speedup computations

public:

	int M;																		///< Number of subbands
	int D;																		///< Downsampling factor
	sbmodel sbm;																///< Subband model

	/// Constructor
	sbconfig( int M_=0, int D_=0 ) : M(M_), D(D_), sbm(M_) {}

	/// Virtual destructor
	virtual
	~sbconfig() {}

	/// Set the analysis filterbank (polyphase)
	virtual
	void set_analysis_fb( const mat<cseq>& H_ )
	{
		assert( (int)H_.size1() == M && (int)H_.size2() == D );
		H = H_;
	}

	/// Set the analysis filterbank (filterbank)
	void set_analysis_fb( const vec<cseq>& h )
	{
		set_analysis_fb( fb2pp_filterbank( h, D ) );
	}

	/// Set the analysis filterbank (complex prototype)
	void set_analysis_fb( const cseq& h0_ )
	{
		_h0         = h0_;
		vec<cseq> h = filterbank_dtft( _h0, M );
		set_analysis_fb( h );
	}

	/// Set the analysis filterbank (real prototype)
	void set_analysis_fb( const rseq& h0 )
	{
		set_analysis_fb( cseq(h0) );
	}

	/// Set the synthesis filterbank (polyphase)
	virtual
	void set_synthesis_fb( const mat<cseq>& F_ )
	{
		assert( (int)F_.size1() == M && (int)F_.size2() == D );
		F  = F_;
		FA = adjoint(F);
	}

	/// Set the synthesis filterbank (filterbank)
	void set_synthesis_fb( const vec<cseq>& f )
	{
		set_synthesis_fb( fb2pp_filterbank( f, D ) );
	}

	/// Set the synthesis filterbank (prototype)
	void set_synthesis_fb( const cseq& f0_ )
	{
		_f0         = f0_;
		vec<cseq> f = filterbank_dtft( _f0, M );
		set_synthesis_fb( f );
	}

	/// Set the synthesis filterbank (real prototype)
	void set_synthesis_fb( const rseq& f0 )
	{
		set_synthesis_fb( cseq(f0) );
	}
	/// Obtain the analysis prototype
	const cseq& h0() const { return _h0; }

	/// Obtain the synthesis prototype
	const cseq& f0() const { return _f0; }

	/// Obtain the polyphase representation of the whole configuration
	mat<cseq> polyphase() const
	{
		mat<cseq> X  = sbm.impulse_response();
		return FA * X * H;
	}

	/// Obtain the cyclic impulse responses
	vec<cseq> impulse_responses() const
	{
		mat<cseq> Gh = polyphase();
		return pp2fb_system( Gh );
	}
};

/// @}
}
}

/// Specialization of matvar for sbindex
template<>
struct matvar<signal::sbsysapp::sbindex> {

	static
	void* create( const string& varname, const signal::sbsysapp::sbindex& x )
	{
		vec<void*> fields(5);
		fields(0) = matvar<double>::create( "m",      x.m );
		fields(1) = matvar<double>::create( "n",      x.n );
		fields(2) = matvar<double>::create( "t",      x.t );
		fields(3) = matvar<double>::create( "real_f", x.real_f );
		fields(4) = matvar<double>::create( "num_f",  x.num_f );
		return matfile::create_struct( varname, fields );
	}

	static
	signal::sbsysapp::sbindex parse( void* pvar )
	{
		void* pm      = matfile::parse_struct( pvar, "m" );
		void* pn      = matfile::parse_struct( pvar, "n" );
		void* pt      = matfile::parse_struct( pvar, "t" );
		void* preal_f = matfile::parse_struct( pvar, "real_f" );
		void* pnum_f  = matfile::parse_struct( pvar, "num_f" );
		signal::sbsysapp::sbindex x;
		x.m      = matvar<double>::parse( pm );
		x.n      = matvar<double>::parse( pn );
		x.t      = matvar<double>::parse( pt );
		x.real_f = matvar<double>::parse( preal_f );
		x.num_f  = matvar<double>::parse( pnum_f );
		return x;
	}
};

/// Specialization of matvar for sbmodel
template<>
struct matvar<signal::sbsysapp::sbmodel> {

	typedef vec<signal::sbsysapp::sbindex> sbivec;

	static
	void* create( const string& varname, const signal::sbsysapp::sbmodel& x )
	{
		vec<void*> fields(4);
		fields(0) = matvar<double>::create( "real_target_f", x.real_target_f );
		fields(1) = matvar<double>::create( "M", x.M );
		fields(2) = matvar<sbivec>::create( "sbindexes", x.sbindexes );
		fields(3) = matvar<rvec>::create( "sbvalues", x.sbvalues );
		return matfile::create_struct( varname, fields );
	}

	static
	signal::sbsysapp::sbmodel parse( void* pvar )
	{
		void* prtf = matfile::parse_struct( pvar, "real_target_f" );
		void* pM   = matfile::parse_struct( pvar, "M" );
		void* psi  = matfile::parse_struct( pvar, "sbindexes" );
		void* psv  = matfile::parse_struct( pvar, "sbvalues" );
		signal::sbsysapp::sbmodel x;
		x.real_target_f = matvar<double>::parse( prtf );
		x.M             = matvar<double>::parse( pM );
		x.sbindexes     = matvar<sbivec>::parse( psi );
		x.sbvalues      = matvar<rvec>::parse( psv );
		return x;
	}
};

}
#endif
