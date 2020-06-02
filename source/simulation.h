#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <limits>
#include <vector>
#include <chrono>

#include <eigen3/Eigen/Core>

// minimal step sizes (except for step=0)
const double MIN_STEP_TIME = 1.0e-9; // step time > 1 ns seems a reasonable minimal time step (as in: no time step will ever be smaller than this)
const double MIN_STEP_POT = 1.0e-6; // step pot > 1 uV seems reasonable (as in: no potential step will ever be smaller than this)

#include "electrodes.h"
#include "system.h"
#include "environment.h"
#include "experiment.h"
#include "simulationcore.h"

class Simulation
{
private:
    std::ostream &output_stream;
    Core core; // simulation core
    Sizing sz; // contains all system size parameters
    double deltaE; // potential step
public:
    double ipc, ipa, Epc, Epa; // for electroanalytical purposes (set during run())

    Environment env;
    Electrode el;
    Experiment exper;
    System sys;
    
    Simulation(std::ostream &output) : output_stream( output ) {}
    ~Simulation() {}

    void setGridSizing(double gamma, double minF, double maxF, double minlograte, double maxlograte)
        { sz.paramGamma = gamma; sz.minF = minF; sz.maxF = maxF; sz.minLogRate = minlograte; sz.maxLogRate = maxlograte; }
    void setPotentialSizing(double deltatheta) { sz.deltaTheta = deltatheta; }
    void setDifferentialOrders(std::size_t numcurr, std::size_t numderiv) { sz.numCurrentPoints = numcurr; sz.numDerivPoints = numderiv; }

    std::size_t run(std::vector<double>&, std::vector<double>&);
private:
    void scanSegment(double, double, bool, std::vector<double>&, std::vector<double>&);
    void delaySegment(double, double);
};

#endif // SIMULATION_H
