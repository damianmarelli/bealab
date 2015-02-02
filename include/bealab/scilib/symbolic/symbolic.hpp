///@file bealab/scilib/symbolic/symbolic.hpp
/// Definitions and methods for symbolic computations.

#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOSYMBOLIC

#ifndef _BEALAB_SYMBOLIC_SYMBOLIC_
#define	_BEALAB_SYMBOLIC_SYMBOLIC_

#include <bealab/scilib/symbolic/maxima.hpp>

namespace bealab
{
namespace symbolic
{
///@defgroup symbolic_symbolic Definitions and methods
/// Definitions and methods for symbolic computations.
///@{

// Import console objects
using bealab::cin;
using bealab::cout;
using bealab::endl;
using bealab::tab;

/// Non-commutative symbol
class symbol_nc : public basic {

	GINAC_DECLARE_REGISTERED_CLASS( symbol_nc, basic );

	string name;																///< Symbol name
	string latex_name;															///< Symbol name in latex

protected:

	/// Print function
	void print( const print_dflt &c, unsigned level = 0 ) const;

	/// Latex print function
	void print_latex( const print_latex &c, unsigned level = 0 ) const;

public:

	/// Constructor
	symbol_nc( const string& n, const string& ln="" ) :
		name(n)
	{
		latex_name = ( ln == "" ) ? n : ln;
	}

	/// Defines the non-comutative nature
	unsigned return_type() const { return return_types::noncommutative; }
};

/// Base class for defining symbolic functions
class symfun_base {

protected:

	unsigned fun_id;															///< Function ID within GiNaC

public:

	/// Function call operator
	template<class... X>
	const function operator()( const X&... x )
	{
		return function( fun_id, ex(x)... );
	}
};

/// Generic symbolic function
class symfun : public symfun_base {
public:

	/// Constructor
	symfun( const std::string& funname, int npar = 1 )
	{
		function_options opts( funname.data(), npar );
		opts.latex_name( funname );
		fun_id = function::register_new( opts );
	}
};

/// Generic non-commutative symbolic function
class symfun_nc : public symfun_base {
public:

	/// Constructor
	symfun_nc( const std::string& funname, int npar = 1 )
	{
		// Return type
		return_type_t return_type = make_return_type_t<symbol_nc>();

		// Set the options
		function_options opts( funname.data(), npar );
		opts.latex_name( funname );
		opts.set_return_type( return_types::noncommutative,
				&return_type );

		// Register the function within GiNaC
		fun_id = function::register_new( opts );
	}
};

/// Show an expression in a PDF document
class latexpdf : protected ofstream {

	string filename;															///< Filename for saving the displayed formulas (and latex code)

public:

	using ofstream::operator<<;

	bool generate_latex = false;												///< Flag to generate latex code.
	bool trace = false;															///< Flag to trace latex output.
	bool interactive = false;													///< If true, an evince instance is launched after destruction.

	/// Constructor
	latexpdf( string fn="bealab" ) :
		ofstream( (fn+".tex").data(), std::ios_base::out ),
		filename(fn)
	{
		auto pbase = dynamic_cast<ofstream*>(this);
		*pbase << latex;
		*pbase << "\\documentclass{article}" << endl;
//		*pbase << "\\pagestyle{empty}" << endl;
		*pbase << "\\usepackage{breqn}" << endl;
		*pbase << "\\usepackage[a4paper]{geometry}" << endl;
//		*pbase << "\\usepackage{color}" << endl;
		*pbase << "\\geometry{tmargin=2cm,bmargin=2cm,lmargin=2cm,rmargin=2cm}" << endl;
		*pbase << "\\begin{document}" << endl;
	}

	/// Destructor
	~latexpdf()
	{
		auto pbase = dynamic_cast<ofstream*>(this);
		*pbase << "\\end{document}" << endl;
		this->close();
		if( trace ) {
			system( "latex " + filename );
			system( "dvipdfm " + filename + " 2>&1" );
		}
		else {
			system( "latex " + filename + " > /dev/null" );
			system( "dvipdfm " + filename + " 2> /dev/null" );
		}
		system( "rm " + filename + ".aux " + filename +".log " + filename +".dvi");
		if( !generate_latex )
			system( "rm " + filename + ".tex" );
		if( interactive )
			system( "evince " + filename + ".pdf 2> /dev/null" );
	}

	/// Handle std::endl
	latexpdf& operator<<( std::ostream& (*e)(std::ostream&) )
	{
		*dynamic_cast<ofstream*>(this) << "\n";
		e(*this);
		return *this;
	}

	/// Output text
	latexpdf& operator<<( const char* ptext )
	{
		*dynamic_cast<ofstream*>(this) << ptext;
		return *this;
	}

//	/// Output text
//	latexpdf& operator<<( const string& text )
//	{
//		*dynamic_cast<ofstream*>(this) << text;
//		return *this;
//	}

	/// Output an equation
	latexpdf& operator<<( const ex& e )
	{
		*dynamic_cast<ofstream*>(this) << "$"<< e << "$";
		return *this;
	}

	/// Output a displayed equation
	latexpdf& equation( const ex& e, const string& varname="" )
	{
		auto pbase = dynamic_cast<ofstream*>(this);
		*pbase << "\\begin{dmath*}" << endl;
		if( varname.size() != 0 )
			*pbase << varname << " = ";
		*pbase << e << endl;
		*pbase << "\\end{dmath*}" << endl;
		return *this;
	}
};

///@name Numeric evaluation

/// Numeric evaluation of an expression
template<class T=double> T numeval( const ex& x );

/// Build a functor<double(double)> via run-time compilation for numeric evaluation of an expression
bealab::function<double(double)> makefun( const ex& fun, const symbol& x );
///@}

///@name Other

/// Create a matrix with non-commutative symbols
matrix symbolic_matrix_nc( int I, int J, const string& var );
///@}

///@name Symbolic functions

/// Error function
const function erf( const ex&e );

/// Matrix transpose
const function trans( const ex&e );
///@}

///@name Calculus

/// Indefinite integral
ex integrate( const ex& f, const ex& x );

/// Definite integral
ex integrate( const ex& f, const ex& x, const ex& a, const ex& b );

/// Limit
ex limit( const ex& f, const ex& x, const ex& val );
///@}

///@name Harmonic analysis

/// Forward Fourier transform
ex fourier_transform( const ex& fun, const realsymbol& t, const realsymbol& f );

/// Inverse Fourier transform
ex ifourier_transform( const ex& fun, const realsymbol& t, const realsymbol& f );

/// Convolution
ex convolution( const ex& f1, const ex& f2, const ex& x );
///@}

///@}
}
}
#endif
#endif
