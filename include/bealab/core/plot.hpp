/// @file bealab/core/plot.hpp
/// Plotting using Gnuplot.

#ifndef _BEALAB_PLOT_
#define	_BEALAB_PLOT_

#include <bealab/scilib/sequences.hpp>
#include <fstream>

namespace bealab
{
/// @defgroup plot Plot
/// Plotting using Gnuplot.
/// @{

/// @name Plotting style attributes

/// One-word attributes (line type, point type and color).
enum plotstyle {
	full, noline, dashed, dotted, dashdotted,						// Line types
	dots, plus, times, star, box, box_filled, circle, circle_filled,
	pyramid, pyramid_filled, ipyramid, ipyramid_filled, diamond,
	diamond_filled, 												// Point types
	black, red, blue, green, magenta, cyan, yellow, grey		 	// Colors
};

/// Attribute for linewidth.
struct linewidth {
	double val;
	linewidth( double v ) : val(v) {}
};

/// Attribute for pointsize.
struct pointsize {
	double val;
	pointsize( double v ) : val(v) {}
};

/// Attribute for pointinterval.
struct pointinterval {
	int val;
	pointinterval( int v ) : val(v) {}
};
/// @}

//------------------------------------------------------------------------------
/// This class models a single curve with its style
class _curve {

	/// @name Curve attributes
	string fname;		///< Filename with the data to plot
	string title;
	bool lines;
	bool points;
	int linetype;
	int pointtype;
	double line_width;
	double point_size;
	int point_interval;
	string color;
	/// @}

	/// @name Parsing of template arguments
	void _parse_args() {}
	template<class... T>
	void _parse_args( const string &str, const T&... sty )
	{
		title = str;
		_parse_args( sty... );
	}
	template<class... T>
	void _parse_args( const linewidth &lw, const T&... sty )
	{
		line_width = lw.val;
		_parse_args( sty... );
	}
	template<class... T>
	void _parse_args( const pointsize &ps, const T&... sty )
	{
		point_size = ps.val;
		_parse_args( sty... );
	}
	template<class... T>
	void _parse_args( const pointinterval &ps, const T&... sty )
	{
		point_interval = ps.val;
		_parse_args( sty... );
	}
	template<class... T>
	void _parse_args( plotstyle ps, const T&... sty )
	{
		switch( ps ) {

		// Linetypes
		case noline:
			lines = false;
			_parse_args( sty... );
			break;
		case full:
			lines    = true;
			linetype = 1;
			_parse_args( sty... );
			break;
		case dashed:
			lines    = true;
			linetype = 2;
			_parse_args( sty... );
			break;
		case dotted:
			lines    = true;
			linetype = 3;
			_parse_args( sty... );
			break;
		case dashdotted:
			lines    = true;
			linetype = 6;
			_parse_args( sty... );
			break;

		// Point types
		case dots:
			points    = true;
			pointtype = 0;
			_parse_args( sty... );
			break;
		case plus:
			points    = true;
			pointtype = 1;
			_parse_args( sty... );
			break;
		case times:
			points    = true;
			pointtype = 2;
			_parse_args( sty... );
			break;
		case star:
			points    = true;
			pointtype = 3;
			_parse_args( sty... );
			break;
		case box:
			points    = true;
			pointtype = 4;
			_parse_args( sty... );
			break;
		case box_filled:
			points    = true;
			pointtype = 5;
			_parse_args( sty... );
			break;
		case circle:
			points    = true;
			pointtype = 6;
			_parse_args( sty... );
			break;
		case circle_filled:
			points    = true;
			pointtype = 7;
			_parse_args( sty... );
			break;
		case pyramid:
			points    = true;
			pointtype = 8;
			_parse_args( sty... );
			break;
		case pyramid_filled:
			points    = true;
			pointtype = 9;
			_parse_args( sty... );
			break;
		case ipyramid:
			points    = true;
			pointtype = 10;
			_parse_args( sty... );
			break;
		case ipyramid_filled:
			points    = true;
			pointtype = 11;
			_parse_args( sty... );
			break;
		case diamond:
			points    = true;
			pointtype = 12;
			_parse_args( sty... );
			break;
		case diamond_filled:
			points    = true;
			pointtype = 13;
			_parse_args( sty... );
			break;

		// Colors
		case black:
			color = "black";
			_parse_args( sty... );
			break;
		case red:
			color = "red";
			_parse_args( sty... );
			break;
		case blue:
			color = "blue";
			_parse_args( sty... );
			break;
		case green:
			color = "green";
			_parse_args(  sty... );
			break;
		case magenta:
			color = "magenta";
			_parse_args( sty... );
			break;
		case cyan:
			color = "cyan";
			_parse_args( sty... );
			break;
		case yellow:
			color = "yellow";
			_parse_args( sty... );
			break;
		case grey:
			color = "grey";
			_parse_args( sty... );
			break;

		// Default
		default:
			_parse_args( sty... );
			break;
		}
	}
	/// @}

	/// Create a temporary filename to hold the plot data
	string create_filename()
	{
		char fn[] = "/tmp/gnuplotXXXXXX";
		int fd = mkstemp( fn );
		if( fd == -1 )
			error("Cannot create temporary file");
		if( close(fd) == -1 )
			error("Cannot close temporary file");
		return fn;
	}

	/// Converts Nan's into Inf's to avoid that gnuplot does crazy things.
	double c2g( double x )
	{
		if( isnan(x) )
			return inf;
		else
			return x;
	}

public:

	/// Constructor for an XY plot (using a variadic template argument)
	template<class... T>
	_curve( const rvec& x, const rvec& y, const T&... sty ) : lines(true),
		points(false), linetype(1), pointtype(1), line_width(1), point_size(1), point_interval(0)
	{
		// Parse template arguments
		_parse_args( sty... );

		// Create temporary filename to hold the plot data
		fname = create_filename();

		// Fill temporary file
		std::ofstream tmp( fname, std::ios_base::out );
		int I = x.size();
		if( I != (int)y.size() )
			error("_curve::_curve() - The x/y-axis vectors do not have the same size");
	    for( int i = 0; i < I; i++ )
	    	tmp << c2g(x(i)) << " " << c2g(y(i)) << endl;
	    tmp.flush();
	    tmp.close();
	}

//	/// Constructor for an XYZ plot (using a variadic template argument)
//	template<class... T>
//	_curve( const rvec& x, const rvec& y, const rvec& z, const T&... sty ) : lines(true),
//		points(false), linetype(1), pointtype(1), line_width(1), point_size(1), point_interval(0)
//	{
//		// Parse template arguments
//		_parse_args( sty... );
//
//		// Create temporary filename to hold the plot data
//		fname = create_filename();
//
//		// Fill temporary file
//		std::ofstream tmp( fname, std::ios_base::out );
//		int I = x.size();
//		if( I != (int)y.size() || I != (int)z.size() )
//			error("_curve::_curve() - The x/y/z-axis vectors do not have the same size");
//	    for( int i = 0; i < I; i++ )
//	    	tmp << c2g(x(i)) << " " << c2g(y(i)) << " " << c2g(z(i)) << endl;
//	    tmp.flush();
//	    tmp.close();
//	}

	/// Constructor for a multi-XY plot (using a variadic template argument)
	template<class... T>
	_curve( const vec<rvec>& x, const vec<rvec>& y, const T&... sty ) : lines(true),
		points(false), linetype(1), pointtype(1), line_width(1), point_size(1), point_interval(0)
	{
		// Parse template arguments
		_parse_args( sty... );

		// Create temporary filename to hold the plot data
		fname = create_filename();

		// Fill temporary file
		std::ofstream tmp( fname, std::ios_base::out );
		int N = x.size();
		if( N != (int)y.size() )
			error("_curve::_curve() - The x/y-axis vectors do not have the same size");
	    for( int n = 0; n < N; n++ ) {
			int I = x(n).size();
			if( I != (int)y(n).size() )
				error("_curve::_curve() - The x/y-axis vectors do not have the same size");
			for( int i = 0; i < I; i++ )
				tmp << c2g(x(n)(i)) << " " << c2g(y(n)(i)) << endl;
			tmp << endl;
			tmp << endl;
	    }
	    tmp.flush();
	    tmp.close();
	}

	/// Constructor for a multi-XYZ plot (using a variadic template argument)
	template<class... T>
	_curve( const vec<rvec>& x, const vec<rvec>& y, const vec<rvec>& z, const T&... sty ) : lines(true),
		points(false), linetype(1), pointtype(1), line_width(1), point_size(1), point_interval(0)
	{
		// Parse template arguments
		_parse_args( sty... );

		// Create temporary filename to hold the plot data
		fname = create_filename();

		// Fill temporary file
		std::ofstream tmp( fname, std::ios_base::out );
		int N = x.size();
		if( N != (int)y.size() || N != (int)z.size() )
			error("_curve::_curve() - The x/y/z-axis vectors do not have the same size");
	    for( int n = 0; n < N; n++ ) {
			int I = x(n).size();
			if( I != (int)y(n).size() || I != (int)z(n).size() )
				error("_curve::_curve() - The x/y/z-axis vectors do not have the same size");
			for( int i = 0; i < I; i++ )
				tmp << c2g(x(n)(i)) << " " << c2g(y(n)(i)) << " " << c2g(z(n)(i)) << endl;
			tmp << endl;
			tmp << endl;
	    }
	    tmp.flush();
	    tmp.close();
	}

	/// Constructor for a matrix plot (using a variadic template argument)
	template<class... T>
	_curve( const rmat& x, const T&... sty ) : lines(false),
		points(false), linetype(1), pointtype(1), line_width(1), point_size(1), point_interval(0)
	{
		// Parse template arguments
		_parse_args( sty... );

		// Create temporary filename to hold the plot data
		fname = create_filename();

		// Fill temporary file
		std::ofstream tmp( fname, std::ios_base::out );
		int I = x.size1();
		int J = x.size2();
		for( int j = 0; j < J; j++ ) {
			for( int i = 0; i < I; i++ )
		    	tmp << j << " " << i << " " << c2g(x(i,j)) << endl;
	    	tmp << endl;
	    }
	    tmp.flush();
	    tmp.close();
	}

	_curve( const _curve& c )
	{
		title          = c.title;
		lines          = c.lines;
		points         = c.points;
		linetype       = c.linetype;
		pointtype      = c.pointtype;
		line_width     = c.line_width;
		point_size     = c.point_size;
		point_interval = c.point_interval;
		color          = c.color;

		fname = create_filename();
		std::ifstream  src(c.fname, std::ios::binary);
		std::ofstream  dst(fname,   std::ios::binary);
		dst << src.rdbuf();
	}

	/// Function called when a figure is cleared
	void clear()
	{
		if( std::remove( fname.data() ) != 0 )
			warning( "Cannot remove temporary file: " + fname );
	}

	/// Get attribute string for the GNUPlot command
	string get_string()
	{
		// Form command line
		ostringstream cmd;
		cmd	<< "'" << fname << "' ";
		cmd << "title '" << title << "' ";
		if( lines && points )
			cmd << "with linespoints ";
		else if( points )
			cmd << "with points ";
		else if( lines )
			cmd << "with lines ";

		if( lines ) {
			cmd << "linetype "  << linetype  << " ";
			cmd << "linewidth " << line_width << " ";
		}
		if( points ) {
			cmd << "pointtype " << pointtype << " ";
			cmd << "pointsize " << point_size << " ";
			cmd << "pointinterval " << point_interval << " ";
		}
		if( color.size() != 0 )
			cmd << "linecolor rgb '" << color << "' ";

		return cmd.str();
	}
};

//------------------------------------------------------------------------------
/// Models one frame with its position and size. Includes a container of _curves
class _frame {
protected:

public:

	enum frmode { nomode, xyplot, xyzplot, surf, map, contour };
	frmode mode;

	/// Curves within this frame
	deque<_curve> curves;

	/// @name Frame properties
	// XXX add Z-properties
	double x_origin, y_origin, x_size, y_size;
	double x_range_1, x_range_2, y_range_1, y_range_2;
	string title, x_label, y_label;
	bool x_log, y_log;
	double x_logbase, y_logbase;
	/// @}

	/// Constructor
	_frame( frmode t=nomode ) : mode(t), x_origin(0), y_origin(0), x_size(1), y_size(1),
			  x_range_1(nan), x_range_2(nan), y_range_1(nan), y_range_2(nan),
			  x_log(false), y_log(false), x_logbase(10), y_logbase(10) {}

	/// Function called when a figure is cleared
	void clear()
	{
		int N = curves.size();
		for( int n = 0; n < N; n++ )
			curves[n].clear();
		curves = deque<_curve>();
	}

	/// Get a command string for the whole frame
	string get_string()
	{
		// Stream used to form commands
		ostringstream cmd;

		// Origin and size
		cmd << "set origin " << x_origin << "," << y_origin << "; ";
		cmd << "set size " << x_size << "," << y_size << "; ";

		// Title
		if( !title.empty() )
			cmd << "set title '" << title << "'; ";

		// Xlabel
		if( !x_label.empty() )
			cmd << "set xlabel '" << x_label << "'; ";

		// Ylabel
		if( !y_label.empty() )
			cmd << "set ylabel '" << y_label << "'; ";

		// Xrange
		if( !isnan(x_range_1) || !isnan(x_range_2) ) {
			cmd << "set xrange [";
			if( !isnan(x_range_1) )
				cmd << x_range_1;
			cmd << ":";
			if( !isnan(x_range_2) )
				cmd << x_range_2;
			cmd <<"]; ";
		}

		// Yrange
		if( !isnan(y_range_1) || !isnan(y_range_2) ) {
			cmd << "set yrange [";
			if( !isnan(y_range_1) )
				cmd << y_range_1;
			cmd << ":";
			if( !isnan(y_range_2) )
				cmd << y_range_2;
			cmd <<"]; ";
		}

		// Logarithmic scaling
		if( x_log )
			cmd << "set logscale x " << x_logbase << ";";
		if( y_log )
			cmd << "set logscale y " << y_logbase << ";";

		// Plot curves
		int N = curves.size();
		for( int n = 0; n < N; n++ ) {
			if( n == 0 )

				switch( mode ) {
				case xyplot:
					cmd << "plot ";
					break;
				case xyzplot:
					cmd << "splot ";
					break;
				case surf:
					cmd << "unset surface; ";
					cmd << "set yrange [:] reverse; ";
					cmd << "set pm3d; ";
					cmd << "splot ";
					break;
				case map:
					cmd << "set view map; ";
					cmd << "set yrange [:] reverse; ";
					cmd << "unset surface; ";
					cmd << "set pm3d; ";
					cmd << "splot ";
					break;
				case contour:
					cmd << "set view map; ";
					cmd << "set yrange [:] reverse; ";
					cmd << "unset surface; ";
					cmd << "unset pm3d; ";
					cmd << "set contour; ";
					cmd << "splot ";
					break;
				default:
					error("_frame mode not set");
					break;
				}

			else
				cmd << "replot ";
			cmd << curves[n].get_string() << "; ";
		}
		if( N > 0 )
			cmd << "clear; replot";

		return cmd.str();
	}
};

//------------------------------------------------------------------------------
/// Low-level driver for a figure, Includes a container of frames.
class figure {

	/// Mapping of terminal names
	std::map<string,string> termname = {
		{"null", "null"},
		{"x11", "x11 noraise"},
		{"qt", "qt noraise"},
		{"eps", "postscript eps enhanced"},
		{"pdf", "pdfcairo"},
		{"svg", "svg"},
		{"png", "png"},
		{"latex", "latex"}
	};

	/// @name Figure attributed
	FILE *pgnuplot;																///< Pointer to a gnuplot session
	deque<_frame> frames;														///< Container of frames within the figure
	string window_title;														///< Window title of filename
	int subplots_size1;															///< Number of rows for an array of subplots
	int subplots_size2;															///< Number of columns for an array of subplots
	bool overlap_f;																///< Overlap flag (if true, all entries of a vector/marix are plotted in the same frame )
//	_frame::frmode spmode;														///< Remembers the type of plot (e.g., xy, xyz) of all frames
	enum terminal_t { t_null, t_x11, t_qt } terminal;
	/// @}

	/// @name Attributes shared by all the frames
	bool x_log;																	///< Flag for logarithmic x-scaling
	bool y_log;																	///< Flag for logarithmic y-scaling
	double x_logbase;															///< Base for logarithmic x-scaling
	double y_logbase;															///< Base for logarithmic y-scaling
	double x_range_1, x_range_2, y_range_1, y_range_2;
	string frame_title, x_label, y_label;
	///@}

	/// Access the frame of a subplot
	_frame* _subplot( int i, int j )
	{
		int J = subplots_size2;
		return &frames[i*J+j];
	}

	void _gnuplot_init()
	{
#ifndef BEALAB_NOGNUPLOT
		// Open pipe
		pgnuplot = popen( "gnuplot 2> /dev/null", "w" );
		if( pgnuplot == NULL )
			error("figure::figure() - Cannot open pipe to gnuplot");

		// Initialize the GNUPlot session
		fputs( "set border 3; set tics nomirror\n", pgnuplot );
		fflush( pgnuplot );

		// Set the terminal
//		x11( false );
		qt();
#endif
	}

	/// Configure the subplots array
	void _set_subplots( int I, int J, _frame::frmode mod=_frame::xyplot )
	{
		// If the subplots array has already been initialized
		if( (subplots_size1 * subplots_size2 ) != 0 ) {

			// If the current sizes are different from the intended ones -> error
			if( subplots_size1 != I || subplots_size2 != J )
				error("figure::set_subplots() - The new sizes are incompatible with the current ones");

			// If the current sizes are equal to the intended ones -> do nothing
			else
				return;
		}

//		// Check if the subplot mode is the right one
//		if( spmode != _frame::nomode && spmode != mod )
//			error("figure::set_subplots() - The mode is incompatible with the current one");
//		else
//			spmode = mod;

		// Create the frames covering the whole screen
		subplots_size1 = I;
		subplots_size2 = J;
		frames.resize(I*J);
		for( int i = 0; i < I; i++ )
			for( int j = 0; j < J; j++ ) {
				frames[i*J+j]           = _frame( mod );
				frames[i*J+j].x_origin  = 1./J * j;
				frames[i*J+j].y_origin  = 1./I * (I-i-1);
				frames[i*J+j].x_size    = 1./J;
				frames[i*J+j].y_size    = 1./I;
				frames[i*J+j].x_log     = x_log;
				frames[i*J+j].y_log     = y_log;
				frames[i*J+j].x_logbase = x_logbase;
				frames[i*J+j].y_logbase = y_logbase;
				frames[i*J+j].x_range_1 = x_range_1;
				frames[i*J+j].x_range_2 = x_range_2;
				frames[i*J+j].y_range_1 = y_range_1;
				frames[i*J+j].y_range_2 = y_range_2;
				frames[i*J+j].title     = frame_title;
				frames[i*J+j].x_label   = x_label;
				frames[i*J+j].y_label   = y_label;
			}
	}

	/// Draw the figure.
	void _draw_force()
	{
		// Number of frames
		int N = frames.size();

		// If more than one frame, enter multiplot mode
		if( N > 1 )
			gnuplot("set multiplot");

		// Plot each frame
		for( int n = 0; n < N; n++ )
			gnuplot( frames[n].get_string() );

		// If more than one frame, quit multiplot mode
		if( N > 1 )
			gnuplot("unset multiplot");
	}

	/// Draw the figure if the terminal is not t_null
	void _draw()
	{
		if( terminal != t_null )
			_draw_force();
	}

	/// Clear the figure's frames, but not the gnuplot window
	void _frames_clear()
	{
		int N = frames.size();
		for( int n = 0; n < N; n++ )
			frames[n].clear();
		frames.resize(0);
		subplots_size1 = 0;
		subplots_size2 = 0;
	}

public:

	bool trace = false;															///< Tracing flag

	/// Constructor
	figure( const string& title="BeaLab" ) : pgnuplot(NULL), window_title(title)
	{
		terminal  = t_qt;
		overlap_f = false;
		x_log     = false;
		y_log     = false;
		x_logbase = 10;
		y_logbase = 10;
		x_range_1 = nan;
		x_range_2 = nan;
		y_range_1 = nan;
		y_range_2 = nan;
		frame_title.clear();
		x_label.clear();
		y_label.clear();
		_frames_clear();
	}

	/// Copy constructor
	figure( const figure& f ) : pgnuplot(NULL)
	{
		// Copy internal variables
		frames         = f.frames;
		window_title   = f.window_title;
		subplots_size1 = f.subplots_size1;
		subplots_size2 = f.subplots_size2;
//		spmode         = f.spmode;
		overlap_f      = f.overlap_f;
		x_log          = f.x_log;
		y_log          = f.y_log;
		x_logbase      = f.x_logbase;
		y_logbase      = f.y_logbase;
		x_range_1      = f.x_range_1;
		x_range_2      = f.x_range_2;
		y_range_1      = f.y_range_1;
		y_range_2      = f.y_range_2;
		frame_title    = f.frame_title;
		x_label        = f.x_label;
		y_label        = f.y_label;
		terminal       = f.terminal;
	}

	/// Destructor
	~figure()
	{
		_frames_clear();
		if( pgnuplot != NULL )
			pclose( pgnuplot );
	}

	/// Put a string in GNUPlot
	figure& gnuplot( const string& cmd_ )
	{
#ifndef BEALAB_NOGNUPLOT
		if( pgnuplot == NULL )
			_gnuplot_init();
		string cmd = cmd_ + "\n";
		if(trace)
			cout << cmd.data() << endl;
		fputs( cmd.data(), pgnuplot );
		fflush( pgnuplot );
#endif
		return *this;
	}

	/// Set the overlap flag (to plot members of containers within the same frame)
	figure& overlap() { overlap_f = true; return *this; }

	/// Clear the figure
	figure& clear()
	{
		// Clear the class frames
		_frames_clear();

		// Clear the gnuplot window
		if( terminal != t_null )
			gnuplot( "clear" );

		return *this;
	}

	/// @name Set frame shared attributes

	/// Set x range
	figure& rangex( double a, double b )
	{
		x_range_1 = a;
		x_range_2 = b;
		return *this;
	}

	/// Set y range
	figure& rangey( double a, double b )
	{
		y_range_1 = a;
		y_range_2 = b;
		return *this;
	}

	/// Set x label
	figure& labelx( const string& label )
	{
		x_label = label;
		return *this;
	}

	/// Set y label
	figure& labely( const string& label )
	{
		y_label = label;
		return *this;
	}

	/// Set title
	figure& title( const string& label )
	{
		frame_title = label;
		return *this;
	}

	/// Set logarithmic x scale
	figure& logx( double base = 10 )
	{
		x_log     = true;
		x_logbase = base;
		return *this;
	}

	/// Set logarithmic y scale
	figure& logy( double base = 10 )
	{
		y_log     = true;
		y_logbase = base;
		return *this;
	}
	/// @}

	/// @name Redirect the output

	/// Redirect the output to /dev/null
	figure& null()
	{
		string output = "set terminal svg; set output '/dev/null'";
		gnuplot( output );
		terminal = t_null;
		return *this;
	}

	/// Redirect the output to an x11 terminal
	figure& x11( bool dforce=true )
	{
	 	string output = "set terminal " + termname["x11"] + " title '" + window_title + "'";
		gnuplot( output );
		if(dforce)
			_draw_force();
		terminal = t_x11;
		return *this;
	}

	/// Redirect the output to a qt terminal
	figure& qt()
	{
	 	string output = "set terminal " + termname["qt"] + " title '" + window_title + "'";
		gnuplot( output );
		_draw_force();
		terminal = t_qt;
		return *this;
	}

	/// Macro to save the plot to a file
#define _SAVE_FORMAT(format) \
	figure& save_##format( string fname=string(), const string& opts=string() ) \
	{ \
		if( terminal == t_null ) \
			_draw_force(); \
		if( fname.empty() ) fname = window_title; \
		string output = "set terminal " + termname[#format] + " " + opts + \
						"; set output '" + fname + "'; replot"; \
		gnuplot( output ); \
		sleep(1); \
		switch( terminal ) { \
		case t_null: \
			null(); \
			break; \
		case t_x11: \
			x11(false); \
			break; \
		case t_qt: \
			qt(); \
			break; \
		} \
		return *this; \
	}
	_SAVE_FORMAT(eps)															///< Save contents to a EPS file
	_SAVE_FORMAT(pdf)															///< Save contents to a PDF file
	_SAVE_FORMAT(svg)															///< Save contents to a SVG file
	_SAVE_FORMAT(png)															///< Save contents to a PNG file
	_SAVE_FORMAT(latex)															///< Save contents to a LATEX file
	/// @}

	/// @name Plot Vectors

	/// Plots a matrix of vectors providing the x-scale
	template<class... T>
	figure& plot( const rvec& x, const mat<rvec> &Y, const T&... sty )
	{
		int I = Y.size1();
		int J = Y.size2();

		if( overlap_f == true )
			_set_subplots( 1, 1 );
		else
			_set_subplots( I, J );

		for( int i = 0; i < I; i++ )
			for( int j = 0; j < J; j++ ) {
				if( overlap_f == true )
					_subplot(0,0)->curves.push_back( _curve( x, Y(i,j), sty... ) );
				else
					_subplot(i,j)->curves.push_back( _curve( x, Y(i,j), sty... ) );
			}
		_draw();
		return *this;
	}

	/// Plots a vector of vectors providing the x-scale
	template<class... T>
	figure& plot( const rvec &x, const vec<rvec> &y, const T&... sty )
	{
		mat<rvec> Y( y.size(), 1 );
		Y.column(0) = y;
		return plot( x, Y, sty... );
	}

	/// Plots a vector providing the x-scale
	template<class... T>
	figure& plot( const rvec &x, const rvec &y, const T&... sty )
	{
		return plot( x, mat<rvec>{{y}}, sty... );
	}

	/// Plots a matrix of vectors
	template<class... T>
	figure& plot( const mat<rvec> &Y, const T&... sty )
	{
		int I = Y.size1();
		int J = Y.size2();

		if( overlap_f == true )
			_set_subplots( 1, 1 );
		else
			_set_subplots( I, J );

		for( int i = 0; i < I; i++ )
			for( int j = 0; j < J; j++ ) {
				if( Y(i,j).size() == 0 )
					continue;
			    rvec x = vrange( 0, Y(i,j).size() );
				if( overlap_f == true )
					_subplot(0,0)->curves.push_back( _curve( x, Y(i,j), sty... ) );
				else
					_subplot(i,j)->curves.push_back( _curve( x, Y(i,j), sty... ) );
			}
		_draw();
		return *this;
	}

	/// Plots a vector of vectors
	template<class... T>
	figure& plot( const vec<rvec> &y, const T&... sty )
	{
		int I = y.size();
		mat<rvec> M(I,1);
		M.column(0) = y;
		return plot( M, sty... );
	}

	/// Plots a vector
	template<class... T>
	figure& plot( const rvec &y, const T&... sty )
	{
		return plot( mat<rvec>{{y}}, sty... );
	}

	/// Plots a matrix of complex vectors
	template<class... T>
	figure& plot( const mat<cvec> &Y, const T&... sty )
	{
		int I = Y.size1();
		int J = Y.size2();

		if( overlap_f == true )
			_set_subplots( 1, 1 );
		else
			_set_subplots( I, J );

		for( int i = 0; i < I; i++ )
			for( int j = 0; j < J; j++ ) {
				if( Y(i,j).size() == 0 )
					continue;
			    rvec x = real( Y(i,j) );
			    rvec y = imag( Y(i,j) );
				if( overlap_f == true )
					_subplot(0,0)->curves.push_back( _curve( x, y, times, noline, sty... ) );
				else
					_subplot(i,j)->curves.push_back( _curve( x, y, times, noline, sty... ) );
			}
		_draw();
		return *this;
	}

	/// Plots a vector of complex vectors
	template<class... T>
	figure& plot( const vec<cvec> &y, const T&... sty )
	{
		return plot( mat<cvec>(y), sty... );
	}

	/// Plots a complex vector
	template<class... T>
	figure& plot( const cvec &y, const T&... sty )
	{
		return plot( mat<cvec>{{y}}, sty... );
	}

	/// Multi-XY plot
	template<class... T>
	figure& plot( const vec<rvec>& x, const vec<rvec> &y, const T&... sty )
	{
		_set_subplots( 1, 1 );
		_subplot(0,0)->curves.push_back( _curve( x, y, sty... ) );
		_draw();
		return *this;
	}
/// @}

	/// @name Plot sequences

	/// Plot a matrix of sequences
	template<class... T>
	figure& plot( const mat<rseq> &X, const T&... sty )
	{
		int I = X.size1();
		int J = X.size2();

		if( overlap_f == true )
			_set_subplots( 1, 1 );
		else
			_set_subplots( I, J );

		for( int i = 0; i < I; i++ )
			for( int j = 0; j < J; j++ ) {
				if( X(i,j).size() == 0 )
					continue;
				rvec t = vrange( X(i,j).t1(), X(i,j).t2()+1 );
				if( overlap_f == true )
					_subplot(0,0)->curves.push_back(
							_curve( t, X(i,j).buffer(), sty... ) );
				else
					_subplot(i,j)->curves.push_back(
							_curve( t, X(i,j).buffer(), sty... ) );
			}
		_draw();
		return *this;
	}

	/// Plot a vector of sequences
	template<class... T>
	figure& plot( const vec<rseq> &Y, const T&... sty )
	{
		mat<rseq> Z( Y.size(), 1 );
		Z.column(0) = Y;
	    return plot( Z, sty... );
	}

	/// Plot a sequence
	template<class... T>
	figure& plot( const rseq &y, const T&... sty )
	{
	    return plot( mat<rseq>{{y}}, sty... );
	}

	/// @name 3D plots

	/// XYZ plot
	template<class... T>
	figure& plot( const rvec &x, const rvec &y, const rvec &z, const T&... sty )
	{
		_set_subplots( 1, 1, _frame::xyzplot );
		_subplot(0,0)->curves.push_back( _curve( vec<rvec>{x}, vec<rvec>{y}, vec<rvec>{z}, sty... ) );
		_draw();
		return *this;
	}

	/// Multi-XYZ plot
	template<class... T>
	figure& plot( const vec<rvec> &x, const vec<rvec> &y, const vec<rvec> &z, const T&... sty )
	{
		_set_subplots( 1, 1, _frame::xyzplot );
		_subplot(0,0)->curves.push_back( _curve( x, y, z, sty... ) );
		_draw();
		return *this;
	}

	/// Surface plot
	template<class... T>
	figure& surf( const rmat &x, const T&... sty )
	{
		_set_subplots( 1, 1, _frame::surf );
		_subplot(0,0)->x_range_1 = 0;
		_subplot(0,0)->x_range_2 = x.size2()-1;
		_subplot(0,0)->y_range_1 = 0;
		_subplot(0,0)->y_range_2 = x.size1()-1;
		_subplot(0,0)->curves.push_back( _curve( x, sty... ) );
		_draw();
		return *this;
	}

	/// Surface plot
	template<class... T>
	figure& map( const rmat &x, const T&... sty )
	{
		_set_subplots( 1, 1, _frame::map );
		_subplot(0,0)->x_range_1 = 0;
		_subplot(0,0)->x_range_2 = x.size2()-1;
		_subplot(0,0)->y_range_1 = 0;
		_subplot(0,0)->y_range_2 = x.size1()-1;
		_subplot(0,0)->curves.push_back( _curve( x, sty... ) );
		_draw();
		return *this;
	}

	/// Surface plot
	template<class... T>
	figure& contour( const rmat &x, const T&... sty )
	{
		_set_subplots( 1, 1, _frame::contour );
		_subplot(0,0)->x_range_1 = 0;
		_subplot(0,0)->x_range_2 = x.size2()-1;
		_subplot(0,0)->y_range_1 = 0;
		_subplot(0,0)->y_range_2 = x.size1()-1;
		_subplot(0,0)->curves.push_back( _curve( x, sty... ) );
		_draw();
		return *this;
	}
	/// @}

	/// @name Shapes

	/// Plot a dot
	template<class... T>
	figure& dot( const rvec& p, const T&... sty )
	{
		int d = p.size();
		assert( d >= 2 );
		assert( d <= 3 );

		if( d == 2 )
			plot( rvec{p(0),p(0)}, rvec{p(1),p(1)}, sty..., noline );
		else
			plot( rvec{p(0),p(0)}, rvec{p(1),p(1)}, rvec{p(2),p(2)}, sty..., noline );
		return *this;
	}

	/// Plot a line
	template<class... T>
	figure& line( const rvec& p1, const rvec& p2, const T&... sty )
	{
		int d = p1.size();
		assert( d == (int)p2.size() );
		assert( d >= 2 );
		assert( d <= 3 );

		if( d == 2 )
			plot( rvec{p1(0),p2(0)}, rvec{p1(1),p2(1)}, sty... );
		else
			plot( rvec{p1(0),p2(0)}, rvec{p1(1),p2(1)}, rvec{p1(2),p2(2)}, sty... );
		return *this;
	}

	/// Plot a legend
	template<class... T>
	figure& legend( const T&... sty )
	{
		plot( rvec{0,0}, rvec{0,0}, rvec{0,0}, sty... );
		return *this;
	}
	/// @}
};

/// @}
}
#endif
