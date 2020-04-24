#include "../source/simulation.h"
#include <stdlib.h>

int main()
{
    // simulation object:
    Simulation* g_sim = new Simulation(std::cout);
    
    // setup CV experiment:
    g_sim->exper->setVertexPotentials({-1.0, 1.0});
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
    
    // OPTIONAL: set advanced settings:
    //g_sim->setGridSizing(_gamma, _minF, _maxF, _minlograte, _maxlograte);
    //g_sim->setPotentialSizing(_deltatheta);
    //g_sim->setDifferentialOrders(static_cast<size_t>(_numcurr), static_cast<size_t>(_numderiv));
    
    // add chemical species to the system (name, initial concentration [mol/m3], diffusion coefficient [m2/s]):
    g_sim->sys->addSpecies( new Species("ox", 1.0, 1.0e-9) );
    g_sim->sys->addSpecies( new Species("red", 0.0, 1.0e-9) );
    g_sim->sys->addSpecies( new Species("product", 0.0, 1.0e-9) );
    
    // add redox steps (oxidized species, reduced species,
    //                  number of electrons, standard potential [V],
    //                  ke [m/s], alpha [-], unused bool)):
    // names should match the ones used in sys->addSpecies(...)
    Redox* redox = new Redox(g_sim->sys->getSpeciesByName("ox"),
                             g_sim->sys->getSpeciesByName("red"),
                             1, -0.5, 1.0, 0.5, false);
    redox->enabled = true;
    g_sim->sys->addRedox(redox);
    
    // add homogeneous reactions of the form A+B <-> C+D (A, B, C, D, k_f, k_b)
    // names should match the ones used in sys->addSpecies(...)
    Reaction* rxn = new Reaction(g_sim->sys->getSpeciesByName("red"),
                                 g_sim->sys->getSpeciesByName(""),
                                 g_sim->sys->getSpeciesByName("product"),
                                 g_sim->sys->getSpeciesByName(""),
                                 10.0, 0.0); // forward, backward rate constant
    rxn->enabled = true;
    g_sim->sys->addReaction(rxn);

    // set up electrode (type and geometry):
    g_sim->el->setType(0); // 0=disk/1=square/2=rectangle/3=cylinder/4=sphere/5=hemisphere
    g_sim->el->setGeom1(1.0e-3); // main size (radius, width)
    g_sim->el->setGeom2(0.0); // (height of rectangle, length of cylinder)

    // run simulation:
    std::vector<double> current, potential;
    size_t sz = g_sim->run(current, potential);
    
    // voltammogram stored in current/potential vectors

    return 0;
}
