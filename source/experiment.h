#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <vector>
#include <iostream>

using namespace std;

extern const double MIN_STEP_TIME;
extern const double MIN_STEP_POT;

class Experiment
{
public:
    double vertexDelay, equilibrationTime, conditioningTime, scanRate;
    double initialPotential, finalPotential, conditioningPotential;
    int numCycles;
    vector<double> vertexPotentials;
public:
    Experiment()
      : vertexDelay( 0.0 ), equilibrationTime( 0.0 ), conditioningTime( 0.0 ), scanRate( 0.0 ),
        initialPotential( 0.0 ), finalPotential( 0.0 ), conditioningPotential( 0.0 ),
        numCycles( 1 ), vertexPotentials( {} )
    {}
    virtual ~Experiment() {}

    double totalTime();

    void setConditioningTime(double t) { conditioningTime = t; }
    void setConditioningPotential(double p) { conditioningPotential = p; }
    void setEquilibration(double t) { equilibrationTime = t; }

    void setInitialPotential(double p) { initialPotential = p; }
    void setVertexPotentials(vector<double> vp) { vertexPotentials = vp; }
    void setFinalPotential(double p) { finalPotential = p; }

    void setVertexDelay(double t) { vertexDelay = t; }

    void setScanRate(double sr) { scanRate = sr; }
    void setNumCycles(int nc) { numCycles = nc; }
};

#endif // EXPERIMENT_H
