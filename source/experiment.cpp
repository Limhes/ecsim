#include <iostream>
#include <cmath>
#include "experiment.h"

double Experiment::totalTime()
{
    double segmentStart, totalTime = 0.0;
    totalTime += conditioningTime;
    totalTime += equilibrationTime;

    for (unsigned int cycle_num = 0; cycle_num < numCycles; cycle_num++) {
        segmentStart = initialPotential;
        for (auto vp: vertexPotentials) {
            totalTime += std::fabs(segmentStart - vp)/scanRate;
            totalTime += vertexDelay;
            segmentStart = vp;
        }
        totalTime += std::fabs(segmentStart - finalPotential)/scanRate;
    }
    totalTime += static_cast<double>(numCycles-1)*vertexDelay;

    return totalTime;
}
