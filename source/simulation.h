#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <limits>
#include <vector>

#include <eigen3/Eigen/Core>

using namespace std;

// minimal step sizes (except for step=0)
const double MIN_STEP_TIME = 1.0e-9; // step time > 1 ns seems a reasonable minimal time step (as in: no time step will ever be smaller than this)
const double MIN_STEP_POT = 1.0e-6; // step pot > 1 uV seems reasonable (as in: no potential step will ever be smaller than this)

#include "electrodes.h"
#include "system.h"
#include "environment.h"
#include "experiment.h"
#include "simulationcore.h"

uint64_t getPosixClockTime();

class Simulation
{
private:
    ostream &output_stream;
    Core core; // simulation core
    Sizing sz; // contains all system size parameters
    double deltaE; // potential step
public:
    double ipc, ipa, Epc, Epa; // for electroanalytical purposes (set during run())

    Environment env;
    Electrode el;
    Experiment exper;
    System sys;
    
    Simulation(ostream &_output) : output_stream( _output ) {}
    ~Simulation() {}

    void setGridSizing(double _gamma, double _minF, double _maxF, double _minlograte, double _maxlograte)
        { sz.paramGamma = _gamma; sz.minF = _minF; sz.maxF = _maxF; sz.minLogRate = _minlograte; sz.maxLogRate = _maxlograte; }
    void setPotentialSizing(double _deltatheta) { sz.deltaTheta = _deltatheta; }
    void setDifferentialOrders(size_t _numcurr, size_t _numderiv) { sz.numCurrentPoints = _numcurr; sz.numDerivPoints = _numderiv; }

    size_t run(vector<double>&, vector<double>&);
private:
    void scanSegment(double, double, bool, vector<double>&, vector<double>&);
    void delaySegment(double, double);
};

#endif // SIMULATION_H
