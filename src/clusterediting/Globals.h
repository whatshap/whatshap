#ifndef GLOBALS_H
#define GLOBALS_H

#include <limits>
#include <cstdlib>

//#include <lemon/time_measure.h>


namespace ysk {

	/**
	 * The time limit in seconds, can be set globally
	 * Note: This is currently only respected by the ILP, the heuristic doesn't care about it and the reduction rules don't take it into account either
	 */
	extern int time_limit;

	extern double threshold;

//	extern lemon::Timer clk;
	extern int verbosity;
	extern int no_threads;
	extern double eps;

	extern bool isTerminated;

	/**
	 * The different states that can occur during the clustering process
	 */
	enum EdgeType{
		UNDECIDED = 0,//!< UNDECIDED The edge is not yet included or excluded from the solution
		PERMANENT = 1,//!< PERMANENT The edge is part of the solution
		FORBIDDEN = 2 //!< FORBIDDEN The edge is not part of the solution
	};


} // namespace ysk

#endif /* GLOBALS_H */

