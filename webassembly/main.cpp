#include "emscripten.h"
#include "simulation.h"
#include <stdlib.h>

extern "C"
{

// simulation data:
double* g_current;
double* g_potential;

// simulation object:
Simulation* g_sim;
std::vector<double> g_vertex_potentials;

EMSCRIPTEN_KEEPALIVE
double* getcurrent()   { return &g_current[0]; }
EMSCRIPTEN_KEEPALIVE
double* getpotential() { return &g_potential[0]; }

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
int addReaction(char* _A, char* _B, char* _C, char* _D, double _kf, double _kb, bool _enabled)
{
    // this is safe, since getSpeciesByName() returns nullptr when string not found,
    // and nullptr leads to "no species" in Reaction constructor:
    Reaction* rxn = new Reaction(g_sim->sys->getSpeciesByName(std::string(_A)),
                                 g_sim->sys->getSpeciesByName(std::string(_B)),
                                 g_sim->sys->getSpeciesByName(std::string(_C)),
                                 g_sim->sys->getSpeciesByName(std::string(_D)),
                                 _kf, _kb);
    rxn->enabled = _enabled;
    g_sim->sys->addReaction(rxn);
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int setElectrode(int _type, double _geom1, double _geom2)
{
    // set up electrode (type and geometry):
    g_sim->el->setType(_type);
    g_sim->el->setGeom1(_geom1);
    g_sim->el->setGeom2(_geom2);
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int setEnvironment(double _temp)
{
    // set up environment:
    g_sim->env->setTemperature(_temp);

    return 0;
}

EMSCRIPTEN_KEEPALIVE
int setCVExperiment(double _tcond, double _econd, double _tequil, double _einit, double _efinal, double _tvertex, double _sr, int _nc)
{
    g_vertex_potentials.clear();
    
    g_sim->exper->setConditioningTime(_tcond);
    g_sim->exper->setConditioningPotential(_econd);
    g_sim->exper->setEquilibration(_tequil);

    g_sim->exper->setInitialPotential(_einit);
    g_sim->exper->setFinalPotential(_efinal);

    g_sim->exper->setVertexDelay(_tvertex);

    g_sim->exper->setScanRate(_sr);
    g_sim->exper->setNumCycles(_nc);
    
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int addCVVertex(double _vp)
{
    g_vertex_potentials.push_back(_vp);
    return 0;
}

EMSCRIPTEN_KEEPALIVE
int setCVVertices()
{
    g_sim->exper->setVertexPotentials(g_vertex_potentials);
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
size_t dosim()
{
    // run simulation:
    std::vector<double> current, potential;
    size_t sz = g_sim->run(current, potential);
    // allocate memory for simulation data:
    g_current =   static_cast<double*>( realloc(g_current,   sz*sizeof(double)) );
    g_potential = static_cast<double*>( realloc(g_potential, sz*sizeof(double)) );
    // copy simulation data to structures accessible in JS:
    for (size_t i = 0; i < sz; i++)
    {
        g_current[i] = current[i];
        g_potential[i] = potential[i];
    }

    // return length of the current/potential arrays:
    return sz;
}


EMSCRIPTEN_KEEPALIVE
int main()
{
    std::cout << "Module loaded." << std::endl;

    // initialize output arrays to nullptr, so realloc knows they're empty:
    g_current = nullptr;
    g_potential = nullptr;

    // setup simulation:
    clearSim();

    return 0;
}

} // end extern C
