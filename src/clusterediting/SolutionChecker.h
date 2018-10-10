#ifndef SOLUTIONCHECKER_H
#define SOLUTIONCHECKER_H

#include "ClusterEditingSolutionLight.h"
#include "Globals.h"
#include <set>

namespace ysk{
    
class SolutionChecker
{
public:
    static const int MARGIN_OF_ERROR = 1;
    
    static bool verifySolution(const LightCompleteGraph& graph, const ClusterEditingSolutionLight& solutions);
};

}

#endif // SOLUTIONCHECKER_H
