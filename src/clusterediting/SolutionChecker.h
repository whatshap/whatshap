#ifndef SOLUTIONCHECKER_H
#define SOLUTIONCHECKER_H

#include "ClusterEditingSolutionLight.h"
#include "Globals.h"
#include <set>
    
class SolutionChecker
{
public:
    static const int MARGIN_OF_ERROR = 1;
    static bool verifySolution ( DynamicSparseGraph& graph, const ClusterEditingSolutionLight& solution );
};

#endif // SOLUTIONCHECKER_H
