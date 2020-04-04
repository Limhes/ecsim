#include <iostream>
#include <cmath>
#include "experiment.h"

using namespace std;

double Experiment::totalTime()
{
    double segmentStart, totalTime = 0.0;
    totalTime += conditioningTime;
    totalTime += equilibrationTime;

    for (int cycle_num = 0; cycle_num < numCycles; cycle_num++) {
        segmentStart = initialPotential;
        for (auto vp: vertexPotentials) {
            totalTime += abs(segmentStart - vp)/scanRate;
            totalTime += vertexDelay;
            segmentStart = vp;
        }
        totalTime += abs(segmentStart - finalPotential)/scanRate;
    }
    totalTime += (numCycles-1)*vertexDelay;

    return totalTime;
}
