#ifndef SRC_COREALGORITHM_H
#define SRC_COREALGORITHM_H

#include <iostream>
#include <vector>
#include <cmath>
#include <bitset>
#include <set>

#include "ClusterEditingSolutionLight.h"
#include "InducedCostHeuristic.h"
#include "StaticSparseGraph.h"
#include "Globals.h"

namespace ysk {

	class CoreAlgorithm{

		public:

		CoreAlgorithm(
				//LightCompleteGraph& graph
				StaticSparseGraph& graph
		)
		:_graph(graph)
		{};

		ClusterEditingSolutionLight run();

		/**
		 * Attempts a "clean" interrupt of the solving process by stopping CPLEX and setting a kill-flag which is checked throughout the process
		 */
		void cancel();

		private:
			//LightCompleteGraph& _graph;
            StaticSparseGraph& _graph;
	};
}

#endif /* SRC_COREALGORITHM_H */
