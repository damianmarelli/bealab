///@file bealab/symbolic/maxima.hpp
/// Interface to call Maxima functions.

#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOSYMBOLIC

#ifndef _BEALAB_SYMBOLIC_MAXIMA_
#define	_BEALAB_SYMBOLIC_MAXIMA_

#include <bealab/core/prelim.hpp>
#include <ginac/ginac.h>

namespace bealab
{
namespace symbolic
{
///@defgroup symbolic_maxima Maxima interface
/// Interface to call Maxima functions.
///@{

// Import GiNaC
using namespace GiNaC;
using GiNaC::function;

///@name Symbols
/// Symbolic special variables defined as in Maxima.
extern symbol plus;																///< To indicate limit from above.
extern symbol minus;															///< To indicate limit from below.
extern symbol inf;																///< Real positive infinity.
extern symbol minf;																///< Real negative infinity.
extern symbol infinity;															///< Complex infinity.
extern symbol und;																///< Undefined result.
extern symbol ind;																///< Bounded indefinite result.
///@}

namespace maxima
{

// Use symbolic::inf instead of bealab::inf
using symbolic::inf;

/// Functor to call functions in Maxima.
/// to cal the Maxima function 'fun':
/// -# Define
/// \verbatim  maxima::functor fun("fun") \endverbatim
/// -# Call (with an arbitrary number of arguments)
/// \verbatim  fun(x,y,z) \endverbatim
class functor {

	/// @name Temporary names
	static const symbol EULER__;
//	const symbol I__     = symbol("%i");
//	const symbol PI__    = symbol("%pi");
	/// @}

	int write_pipe;																///< Pipe to write to Maxima
	int read_pipe;																///< Pipe to read from Maxima
	string _fname;																///< Function name in maxima
	symtab _table;																///< Table of symbols to parse a maxima expression

	/// Replace all the appearances of 'oldpat' by 'newpat' in 'str'.
	string str_replace( string str, const string& oldpat, const string& newpat )
	{
		int n = 0;
		while(true) {
			n = str.find( oldpat, n );
			if( n == -1 )
				break;
			str.replace( n, oldpat.size(), newpat );
			n += newpat.size();
		};
		return str;
	}

	/// Send a command to maxima
	void maxima_send( const string& command )
	{
		// Put command in maxima
		uint nwritten = write( write_pipe, command.data(), command.size()+1 );
		if( nwritten != command.size()+1 )
			error("maxima::functor::maxima_send_receive() - Writing to Maxima failed");

		// Trace
		if(trace)
			cout << "Command sent to Maxima: " << endl
			     << command << endl << endl;
	}

	/// Read a result from maxima
	string maxima_receive( int nans = INT_MAX )
	{
		// Trace
		if(trace)
			cout << "Answer(s) received from Maxima: " << endl;

		// Read result from maxima
		string str;
		int N = 256;
		char buff[N+1];
		for( int a = 0; a < nans; a++ ) {
			int n = read( read_pipe, buff, N );
			if( n == 0 )
				break;
			buff[n] = 0;
			string sbuff(buff);
			str += sbuff;

			// Trace
			if( trace )
				cout << sbuff;
		}


		// Return the result
		return str;
	}

	string diff2string( const ex& e, const deque<int>& idxs )
	{
		ostringstream rv;
		rv << "diff(";
		if( idxs.size() > 1 ) {
			deque<int> idxsr( idxs.begin(), idxs.end()-1 );
			rv << diff2string( e, idxsr );
		}
		else {
			rv << ex_to<function>(e).get_name() << "(";
			int I = e.nops();
			for( int i = 0; i < I; i++ ) {
				rv << ginac2maxima( e.op(i) );
				if( i < I-1 )
					rv << ",";
			}
			rv << ")";
		}
		rv << "," << ginac2maxima( e.op(idxs.back()) ) << ")";
		return rv.str();
	}

	/// Converts a GiNaC expression into a Maxima string
	string ginac2maxima( const ex& e )
	{
//		cout << ex_to<basic>(e).class_name() << endl;

		// Stream to collect the output
		ostringstream rv;

		// Number of operands
		int I = e.nops();

		// Non compound class
		if( I == 0 ) {
			if( is_a<constant>(e) ) {
				if( e.is_equal(I) )
					rv << "%i";
				else if( e.is_equal(Pi) )
					rv << "%pi";
			}
			else
				rv << e;
		}

		// Derivative
		else if( is_a<fderivative>(e) ) {

			// Parse
			ostringstream os;
			os << e;
		    std::istringstream is(os.str());
		    is.get();
		    is.get();
		    deque<int> idxs;
			while( true ) {
				int idx;
				is >> idx;
				idxs.push_back(idx);
				if( is.get() == ']' )
					break;
			}

			// Output
			rv << diff2string( e, idxs );
//			rv << "diff(" << ex_to<function>(e).get_name() << "(";
//			for( int i = 0; i < I; i++ ) {
//				rv << ginac2maxima( e.op(i) );
//				if( i < I-1 )
//					rv << ",";
//			}
//			rv << ")," << ginac2maxima( e.op(0) ) << ")";
		}

		// Function
		else if( is_a<function>(e) ) {
			rv << ex_to<function>(e).get_name() << "(";
			for( int i = 0; i < I; i++ ) {
				rv << ginac2maxima( e.op(i) );
				if( i < I-1 )
					rv << ",";
			}
			rv << ")";
		}

		// Addition
		else if( is_a<add>(e) ) {
			for( int i = 0; i < I; i++ ) {
				rv << ginac2maxima( e.op(i) );
				if( i < I-1 )
					rv << "+";
			}
		}

		// Multiplication
		else if( is_a<mul>(e) ) {
			for( int i = 0; i < I; i++ ) {
				rv << "(" << ginac2maxima( e.op(i) ) << ")";
				if( i < I-1 )
					rv << "*";
			}
		}

		// Power
		else if( is_a<power>(e) )
			rv << "(" << ginac2maxima( e.op(0) ) << ")^("
			   << ginac2maxima( e.op(1) ) << ")";

		return rv.str();
	}

	/// Converts a Maxima string into a GiNaC expression
	ex maxima2ginac( string maxex )
	{
		// Handle Maxima questions
		if( maxex.find('?') != string::npos )
			error( "Maxima asked: " + maxex );

		// Handle divergent integrals
		if( maxex.find("defint: integral is divergent.") != string::npos )
			error( "Maxima says: integral is divergent" );

		// Remove '
		maxex = str_replace( maxex, "'", "" );

		// Replace %pi -> Pi
		maxex = str_replace( maxex, "%pi", "Pi" );

		// Replace %i -> I
		maxex = str_replace( maxex, "%i", "I" );

		// Replace %e-( by %e(-
		maxex = str_replace( maxex, "%e^-(", "%e^(-" );

		// Replace %e by the temporal symbol EULER__
		maxex = str_replace( maxex, "%e", EULER__.get_name() );

		// Parse the result
		parser reader( _table );
		reader.strict = true;
		ex e = reader( maxex );


		// Remove temporary notation
		e = subs( e, pow(EULER__,wild()) == exp(wild()) );

		return e;
	}

	/// Get all the symbols in an expression
	set<ex,ex_is_less> get_symbols( const ex& e )
	{
		set<ex,ex_is_less> s;
		int I = e.nops();
		if( I )
			for( int i = 0; i < I; i++ ) {
				set<ex,ex_is_less> si = get_symbols( e.op(i) );
				s.insert( si.begin(), si.end() );
			}
		else
			if( is_a<symbol>(e) )
				s.insert( e );
		return s;
	}

	/// Get all the symbols in a vector of expression
	set<ex,ex_is_less> get_symbols( const vector<ex>& exs )
	{
		set<ex,ex_is_less> s;
		int I = exs.size();
		for( int i = 0; i < I; i++ ) {
			set<ex,ex_is_less> si = get_symbols( exs[i] );
			s.insert( si.begin(), si.end() );
		}
		return s;
	}

	/// Fill the table of symbols in the expressions sent as parameters
	void fill_symtab( const vector<ex>& params )
	{
		auto syms = get_symbols( params );
		symtab table;
		for( auto it = syms.begin(); it != syms.end(); it++ )
			table[ex_to<symbol>(*it).get_name()] = *it;
		table[EULER__.get_name()]  = EULER__;
		table[plus.get_name()]     = plus;
		table[minus.get_name()]    = minus;
		table[inf.get_name()]      = inf;
		table[minf.get_name()]     = minf;
		table[und.get_name()]      = und;
		table[ind.get_name()]      = ind;
		table[infinity.get_name()] = infinity;
		_table = table;
	}

public:

	static bool trace;

	/// Constructor
	functor( const string& fname ) : _fname(fname) {}

	/// Call the maxima function
	template<class... A>
	ex operator()( const A&... a )
	{
		// Parse the variadic template arguments
		deque<ex> params_ = parse_variadic_template<ex>( a... );
		vector<ex> params( params_.begin(), params_.end() );

		// Build symbolic table
		fill_symtab( params );

		// Start Maxima
		rw_pipe( "/usr/bin/maxima", {"--very-quiet"}, &write_pipe, &read_pipe );

		// Command to send
		ostringstream cmd;
		int nans = 0;															// Number of answers to expect from maxima

		// Set 2d output
		cmd << "display2d:false$ ";
		cmd << "ratprint:false$ ";

//		// Set conversions
//		cmd << "Pi : %pi; I : %i; ";
//		A += 2;

		// Set assumptions
		for( auto it = _table.begin(); it != _table.end(); it++ )
			if( is_a<possymbol>(it->second) ) {
				cmd << "assume( " << ex_to<symbol>(it->second).get_name()
					<< " > 0 ); ";
				nans++;
			}
			else if( is_a<realsymbol>(it->second) ) {
				cmd << "declare( " << ex_to<symbol>(it->second).get_name()
					<< ", real ); ";
				nans++;
			}

		// Maxima function call
		cmd << _fname << "( ";
		int P = params.size();
		for( int p = 0; p < P-1; p++ )
			cmd << ginac2maxima(params[p]) << ", ";
		cmd << ginac2maxima(params[P-1]) << " ); ";

		// Maxima quit() function
		cmd << "quit();";

		// Send command to Maxima
		maxima_send( cmd.str() );

		// Receive answers from Maxima
		maxima_receive( nans );
		string result = maxima_receive();

		// Parse the result and return
		return maxima2ginac( result );
	}
};

}
///@}
}
}
#endif
#endif

/*
 * TI VOGLIO TANTO BENE :* (smak!)
 */
