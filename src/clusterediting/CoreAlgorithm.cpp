#include "CoreAlgorithm.h"

using namespace std;

/**
    * Generates a set of solutions from its internal parameters (instance and parameter set)
    * @return The solution set as CES*
    */
ClusterEditingSolutionLight CoreAlgorithm::run() {
    //At the beginning of a run we set the global termination flag to false ->
    //If we use this software as a lib it might have been set to true in a previous run
    isTerminated = false;

    //Just a killswitch to prevent the program from running if the user has already cancelled it at this point
    if (isTerminated){
        ClusterEditingSolutionLight sol;
        return sol;
    }
    
    // Run heuristic
    if (verbosity > 1) {
        cout << "Starting CE Heuristic!" << endl;
    }
    
    InducedCostHeuristic hl(_graph, bundleEdges);
    ClusterEditingSolutionLight sol = hl.solve();
    
    if (verbosity > 2) {
        cout << "Total editing cost:\t" << sol.getTotalCost() << endl;
    }

    return sol;
}

void CoreAlgorithm::cancel(){
    isTerminated = true;
}
