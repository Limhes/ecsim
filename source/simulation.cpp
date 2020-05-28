#include <iostream>
#include "simulation.h"

uint64_t getPosixClockTime()
{
    struct timespec ts = {};
    return (clock_gettime (CLOCK_MONOTONIC, &ts) == 0) ? static_cast<uint64_t>(ts.tv_sec * 1e9 + ts.tv_nsec) : 0;
}

/*===============================================================================================
 * HIGHEST LEVEL CLASS METHODS
 *=============================================================================================*/

size_t Simulation::run(vector<double> &current, vector<double> &potential)
{
    // finalize system (variable normalization, thermal equilibration):
    sys.finalize(el.epsilon);
    // initialize system sizing (determine all sizing parameter for the simulation):
    deltaE = sz.initialize(&sys, &el, &env, &exper, output_stream);
    // initialize simulation core ():
    core.initialize(&sz, &sys, exper.initialPotential);

    // initialize electroanalytical parameters:
    ipc = 0.0;
    ipa = 0.0;
    Epc = 0.0;
    Epa = 0.0;

    uint64_t startTime, endTime; // timing variables [ns]
    startTime = getPosixClockTime();
    
    // conditioning step:
    delaySegment(exper.conditioningPotential, exper.conditioningTime);
    
    // equilibration step:
    delaySegment(exper.initialPotential, exper.equilibrationTime);
    
    // scan cycles:
    double segmentStartPotential; // start potential of scan segment
    bool recordCurrent; // record current in segment?
    for (int cycle = 0; cycle < exper.numCycles; cycle++) {
        recordCurrent = (cycle == exper.numCycles - 1); // warning: hard-coded option to record current only of last cycle

        segmentStartPotential = exper.initialPotential;
        for (auto vp: exper.vertexPotentials) {
            // scan up to vertex:
            scanSegment(segmentStartPotential, vp, recordCurrent, current, potential);
            segmentStartPotential = vp;
            // vertex delay:
            delaySegment(vp, exper.vertexDelay);
        }
        // scan from last vertex to final potential:
        scanSegment(segmentStartPotential, exper.finalPotential, recordCurrent, current, potential);
    }
    
    endTime = getPosixClockTime();
    output_stream << "Est. simulation time (" << current.size() << " steps &amp; " << sz.numGridPoints << " grid points) = " << static_cast<double>(endTime - startTime)/1.0e9 << " s<br />" << endl;

    return current.size();
}

void Simulation::delaySegment(double pot, double time)
{
    if (time > MIN_STEP_TIME)
    {
        int num_points = static_cast<int>( time * exper.scanRate / deltaE );
        
        for (int d = 0; d < num_points; d++)
            core.solveSystem(pot);
    }
}

void Simulation::scanSegment(double potential_start, double potential_stop, bool record_current, vector<double> &current, vector<double> &potential)
{
    double curr;

    double sign = 0.0;
    if (potential_start < potential_stop) sign = 1.0;
    else if (potential_start > potential_stop) sign = -1.0;

    int num_points = static_cast<int>(abs(potential_start-potential_stop)/deltaE);
    double pot = potential_start;
    for (int d = 0; d < num_points; d++)
    {
        // do magic:
        core.solveSystem(pot);

        if (record_current) {
            // get current:
            curr = core.calcCurrentFromFlux();

            // update electroanalytical metrics:
            if (curr > ipa) { ipa = curr; Epa = pot; }
            if (curr < ipc) { ipc = curr; Epc = pot; }

            // add data points to current & potential vectors:
            potential.push_back(pot);
            current.push_back(curr);
        }

        // increase potential:
        pot += sign*deltaE;
    }
}
