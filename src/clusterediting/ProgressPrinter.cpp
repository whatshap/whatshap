#include "ProgressPrinter.h"
#include <iostream>
#include <cmath>

ProgressPrinter::ProgressPrinter(std::string name, uint64_t initialProgress, uint64_t totalSteps) :
    name(name),
    progress(initialProgress),
    totalSteps(totalSteps),
    percentage(0),
    progressOfNextPercentage(0) 
{
//     checkProgress();
}

void ProgressPrinter::step() {
    progress++;
    checkProgress();
}

void ProgressPrinter::setProgress(uint64_t newProgress) {
    progress = std::max(progress, std::min(totalSteps, newProgress));
    checkProgress();
}

void ProgressPrinter::setFinished() {
    progress = totalSteps;
    checkProgress();
}

std::string ProgressPrinter::getName() {
    return name;
}

uint64_t ProgressPrinter::getProgress() {
    return progress;
}

uint64_t ProgressPrinter::getTotalSteps() {
    return totalSteps;
}

void ProgressPrinter::checkProgress() {
    if (progress >= progressOfNextPercentage) {
        percentage = totalSteps != 0 ? (progress*100)/totalSteps : 100;
        if (verbosity >= 1) {
//             if (verbosity >= 2 || progress == totalSteps) {
//                 std::cout<<name<<" .. "<<percentage<<"%"<<std::endl;
//             } else {
//                 std::cout<<name<<" .. "<<percentage<<"%\r"<<std::flush;;
//             }
            std::cout<<name<<" .. "<<percentage<<"%\r"<<std::flush;;
        }
        progressOfNextPercentage = (uint64_t)std::ceil(((percentage+1)*totalSteps) / 100.0);
    }
}
