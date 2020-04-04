#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include <iostream>

using namespace std;

const string NO_SPECIES_NAME = "(no species)";

class Species {
public:
    double initialConcentration, normalizedConcentration, normalizedAndEquilibratedConcentration;
    double diffusionConstant, normalizedDiffusionConstant;
    double Del, conc; // for equilibration purposes
    size_t index;
    string name;
    bool inRedox, inReaction;

    Species(string n, double c, double d)
      : initialConcentration( c ), normalizedConcentration( 0.0 ), normalizedAndEquilibratedConcentration( 0.0 ),
        diffusionConstant( d ), normalizedDiffusionConstant( 0.0 ), Del( 0.0 ), conc( 0.0 ),
        index( 0 ), name( n ), inRedox( false ), inReaction( false )
    {}
};

class Reaction
{ // specLHS1 + specLHS2 <==> specRHS1 + specRHS2
public:
    Species *specLHS1, *specLHS2, *specRHS1, *specRHS2;
    double rateConstantForward, rateConstantForwardNormalized; // holds the normalized forward rate constant: Kf * t (unimolecular) or Kf * t * max_conc (bimolecular)
    double rateConstantBackward, rateConstantBackwardNormalized; // holds the normalized backward rate constant: Kb * t (unimolecular) or Kb * t * max_conc (bimolecular)
    size_t index;
    bool enabled;

    Reaction(Species *a, Species *b, Species *c, Species *d, double kf, double kb)
      : specLHS1( a ), specLHS2( b ), specRHS1( c ), specRHS2( d ),
        rateConstantForward( kf ), rateConstantForwardNormalized( 0.0 ),
        rateConstantBackward( kb ), rateConstantBackwardNormalized( 0.0 ), index( 0 ), enabled( false )
    {}
};

class Redox
{ // specOxidized + n*e- <==> specReduced
public:
    Species *specOxidized, *specReduced;
    int numberElectrons; // holds the number of electrons
    double standardPotential; // holds the potential [V]
    double rateConstantHetero, rateConstantHeteroNormalized; // holds the normalized heterogeneous rate constant: Ke * delta_t / delta_x
    double alpha; // holds the alpha coefficient
    bool isRedoxReversible; // if true, this step is Nernstian
    size_t index;
    bool enabled;

    Redox(Species *ox, Species *red, int n, double E, double Ke, double a, bool rev)
      : specOxidized( ox ), specReduced( red ), numberElectrons( n ), standardPotential( E ),
        rateConstantHetero( Ke ), rateConstantHeteroNormalized( 0.0 ), alpha( a ), isRedoxReversible( rev ), index( 0 ), enabled( false )
    {}
};

class System
{
public:
    // reaction network structs (vecX is a subset of vecAllX, where vecX contains the items used in the simulation):
    vector<Species*> vecSpecies, vecAllSpecies; // vector holding a s_species struct per species in the system
    vector<Reaction*> vecReactions, vecAllReactions; // vector holding a s_reaction struct per homogeneous reaction in the system
    vector<Redox*> vecRedox, vecAllRedox; // vector holding a s_redox struct per redox couple in the system
    // system metrics:
    double maxConcentration, maxDiffusionConstant, maxRateConstantChem; // Maximum concentration [mol/m3], diffusion coefficient [m2/s], and homo/hetero max rate constants [1/s] and [??]
public:
    System() : maxConcentration( 0.0 ), maxDiffusionConstant( 0.0 ), maxRateConstantChem( 0.0 ) {}
    ~System();

    vector<double>::size_type addSpecies(Species*);
    vector<double>::size_type addRedox(Redox*);
    vector<double>::size_type addReaction(Reaction*);

    string generateUniqueSpeciesName();
    bool isSpeciesPresentWithName(string);
    Species* getSpeciesByName(string);

    void finalize(double); // parameter: characteristic dimension of system (epsilon)
private:
    void addSpeciesToSystem(Species*);

    void setActiveSystem();
    void updateSpeciesProperties();
    void calcMaxRateConstants(double); // parameter: characteristic dimension of system (epsilon)
    int equilibrateConcentrations();
};

#endif // SYSTEM_H
