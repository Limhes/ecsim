#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <vector>
#include <iostream>
#include <QObject>

using namespace std;

extern const double MIN_STEP_TIME;
extern const double MIN_STEP_POT;

class Experiment : public QObject {
    Q_OBJECT

public:
    double vertexDelay, equilibrationTime, conditioningTime, scanRate;
    double initialPotential, finalPotential, conditioningPotential;
    int numCycles;
    vector<double> vertexPotentials;
public:
    Experiment(double IP, const vector<double>& VP, double FP, double sr)
      : vertexDelay( 0.0 ), equilibrationTime( 0.0 ), conditioningTime( 0.0 ), scanRate( sr ),
        initialPotential( IP ), finalPotential( FP ), conditioningPotential( 0.0 ),
        numCycles( 1 ), vertexPotentials( VP )
    {}
    virtual ~Experiment() {}

    double totalTime();
public slots:
    void setConditioningTime(double t) { conditioningTime = t; }
    void setConditioningPotential(double p) { conditioningPotential = p; }
    void setEquilibration(double t) { equilibrationTime = t; }

    void setInitialPotential(double p) { initialPotential = p; }
    void setFinalPotential(double p) { finalPotential = p; }

    void setVertexDelay(double t) { vertexDelay = t; }

    void setScanRate(double sr) { scanRate = sr; }
    void setNumCycles(int nc) { numCycles = nc; }
};

#endif // EXPERIMENT_H
