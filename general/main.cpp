#include "simulation.h"
#include <stdlib.h>

EMSCRIPTEN_KEEPALIVE
void clearSim()
{
    delete g_sim;
    g_sim = new Simulation(std::cout);

    std::cout << "Simulation object created." << std::endl;
}

EMSCRIPTEN_KEEPALIVE
int addSpecies(char* _name, double _conc, double _diff)
{
    std::cout << _name << " " << _conc << " " << _diff << std::endl;
    g_sim->sys->addSpecies( new Species(std::string(_name), _conc, _diff) );
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int addRedox(char* _ox, char* _red, int _n, double _pot, double _ke, double _alpha, bool _rev, bool _enabled)
{
    Redox* redox = new Redox(g_sim->sys->getSpeciesByName(std::string(_ox)),
                             g_sim->sys->getSpeciesByName(std::string(_red)),
                             _n, _pot, _ke, _alpha, _rev);
    redox->enabled = _enabled;
    g_sim->sys->addRedox(redox);
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int setGridSizing(double _gamma, double _minF, double _maxF, double _minlograte, double _maxlograte)
{
    g_sim->setGridSizing(_gamma, _minF, _maxF, _minlograte, _maxlograte);
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int setPotentialSizing(double _deltatheta)
{
    g_sim->setPotentialSizing(_deltatheta);
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int setDifferentialOrders(int _numcurr, int _numderiv)
{
    g_sim->setDifferentialOrders(static_cast<size_t>(_numcurr), static_cast<size_t>(_numderiv));
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int main()
{
    // simulation object:
    Simulation* g_sim = new Simulation(std::cout);
    
    // setup CV experiment:
    g_sim->exper->setVertexPotentials([-1.0, 1.0]);
    g_sim->exper->setConditioningTime(0.0);
    g_sim->exper->setConditioningPotential(0.0);
    g_sim->exper->setEquilibration(0.0);
    g_sim->exper->setInitialPotential(0.0);
    g_sim->exper->setFinalPotential(0.0);
    g_sim->exper->setVertexDelay(0.0);
    g_sim->exper->setScanRate(0.1);
    g_sim->exper->setNumCycles(1);
    
    // set up environment:
    g_sim->env->setTemperature(293.15); // [K]
    
    
    // this is safe, since getSpeciesByName() returns nullptr when string not found,
    // and nullptr leads to "no species" in Reaction constructor:
    Reaction* rxn = new Reaction(g_sim->sys->getSpeciesByName(std::string('A')),
                                 g_sim->sys->getSpeciesByName(std::string('')),
                                 g_sim->sys->getSpeciesByName(std::string('B')),
                                 g_sim->sys->getSpeciesByName(std::string('')),
                                 _kf, _kb); // forward, backward rate constant
    rxn->enabled = true;
    g_sim->sys->addReaction(rxn);

    // set up electrode (type and geometry):
    g_sim->el->setType(_type);
    g_sim->el->setGeom1(_geom1);
    g_sim->el->setGeom2(_geom2);

    // run simulation:
    std::vector<double> current, potential;
    size_t sz = g_sim->run(current, potential);
    
    

    return 0;
}
