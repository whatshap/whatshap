#ifndef PROGRESSPRINTER_H
#define PROGRESSPRINTER_H

#include <string>
#include "Globals.h"

class ProgressPrinter {

public:
    ProgressPrinter(std::string name, uint64_t initialProgress, uint64_t totalSteps);
    void step();
    void setProgress(uint64_t newProgress);
    void setFinished();
    std::string getName();
    uint64_t getProgress();
    uint64_t getTotalSteps();


private:
    std::string name;
    uint64_t progress;
    uint64_t totalSteps;
    uint64_t percentage;
    uint64_t progressOfNextPercentage;
    void checkProgress();
};

#endif // PROGRESSPRINTER_H
