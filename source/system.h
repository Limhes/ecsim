#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include <iostream>

using namespace std;

const string NO_SPECIES_NAME = "(no species)";

class Species {
private:
    double initialConcentration, normalizedConcentration, normalizedAndEquilibratedConcentration;
    double diffusionConstant, normalizedDiffusionConstant;
    size_t index;
    string name;
public:
    double getConcInit() const { return initialConcentration; }
    void setConcInit(double c) { initialConcentration = c; }
    double getConcEquil() const { return normalizedAndEquilibratedConcentration*initialConcentration/normalizedConcentration; } // TODO: catch divby0
    double getDiff() const { return diffusionConstant; }
    void setDiff(double d) { diffusionConstant = d; }
    
    double getConcNorm() const { return normalizedConcentration; }
    double getConcNormEquil() const { return normalizedAndEquilibratedConcentration; }
    void setConcNormEquil(double c) { normalizedAndEquilibratedConcentration = c; }
    double getDiffNorm() const { return normalizedDiffusionConstant; }
    void normalizeConc(double conc_max) { normalizedConcentration = initialConcentration/conc_max; } // TODO: catch divby0
    void normalizeDiff(double diff_max) { normalizedDiffusionConstant = diffusionConstant/diff_max; } // TODO: catch divby0
    
    void setName(string new_name) { name = new_name; }
    string getName() const { return name; }
    void setIndex(size_t new_index) { index = new_index; }
    size_t getIndex() const { return index; }

    Species(string n, double c, double d)
      : initialConcentration( c ), normalizedConcentration( 0.0 ), normalizedAndEquilibratedConcentration( 0.0 ),
        diffusionConstant( d ), normalizedDiffusionConstant( 0.0 ), index( 0 ), name( n )
    {}
};

class Reaction
{ // specLHS1 + specLHS2 <==> specRHS1 + specRHS2
    Species *specLHS1, *specLHS2, *specRHS1, *specRHS2;
    double rateConstantForward, rateConstantForwardNormalized; // holds the normalized forward rate constant: Kf * t (unimolecular) or Kf * t * max_conc (bimolecular)
    double rateConstantBackward, rateConstantBackwardNormalized; // holds the normalized backward rate constant: Kb * t (unimolecular) or Kb * t * max_conc (bimolecular)
    size_t index;
    bool enabled;
public:
    Species* getSpecLHS1() const { return specLHS1; }
    Species* getSpecLHS2() const { return specLHS2; }
    Species* getSpecRHS1() const { return specRHS1; }
    Species* getSpecRHS2() const { return specRHS2; }
    
    double getKf() const { return rateConstantForward; }
    double getKb() const { return rateConstantBackward; }
    void setKf(double k) { rateConstantForward = k; }
    void setKb(double k) { rateConstantBackward = k; }
    double getKfNorm() const { return rateConstantForwardNormalized; }
    double getKbNorm() const { return rateConstantBackwardNormalized; }
    void normalizeKf(double norm) { rateConstantForwardNormalized = rateConstantForward / norm; } // TODO: catch divby0
    void normalizeKb(double norm) { rateConstantBackwardNormalized = rateConstantBackward / norm; } // TODO: catch divby0
    
    void setIndex(size_t new_index) { index = new_index; }
    size_t getIndex() const { return index; }
    
    bool isEnabled() const { return enabled; }
    void enable() { enabled = true; }
    void disable() { enabled = false; }
    
    Reaction(Species* a, Species* b, Species* c, Species* d, double kf, double kb)
      : specLHS1( a ), specLHS2( b ), specRHS1( c ), specRHS2( d ),
        rateConstantForward( kf ), rateConstantForwardNormalized( 0.0 ),
        rateConstantBackward( kb ), rateConstantBackwardNormalized( 0.0 ), index( 0 ), enabled( false )
    {}
};

class Redox
{ // specOxidized + n*e- <==> specReduced
    Species *specOxidized, *specReduced;
    int numberElectrons; // holds the number of electrons
    double standardPotential; // holds the potential [V]
    double rateConstantHetero, rateConstantHeteroNormalized; // holds the normalized heterogeneous rate constant: Ke * delta_t / delta_x
    double alpha; // holds the alpha coefficient
    size_t index;
    bool enabled;
public:
    Species* getSpecOx() const { return specOxidized; }
    Species* getSpecRed() const { return specReduced; }
    
    double getKe() const { return rateConstantHetero; }
    void setKe(double k) { rateConstantHetero = k; }
    double getKeNorm() const { return rateConstantHeteroNormalized; }
    void normalizeKe(double norm) { rateConstantHeteroNormalized = rateConstantHetero / norm; } // TODO: catch divby0
    double getE0() const { return standardPotential; }
    void setE0(double E) { standardPotential = E; }
    int getNe() const { return numberElectrons; }
    void setNe(int n) { numberElectrons = n; }
    double getAlpha() const { return alpha; }
    void setAlpha(double a) { alpha = a; }
    
    void setIndex(size_t new_index) { index = new_index; }
    size_t getIndex() const { return index; }
    
    bool isEnabled() const { return enabled; }
    void enable() { enabled = true; }
    void disable() { enabled = false; }

    Redox(Species* ox, Species* red, int n, double E, double Ke, double a)
      : specOxidized( ox ), specReduced( red ), numberElectrons( n ), standardPotential( E ),
        rateConstantHetero( Ke ), rateConstantHeteroNormalized( 0.0 ), alpha( a ), index( 0 ), enabled( false )
    {}
};

class System
{
private:
    // system metrics:
    double maxConcentration, maxDiffusionConstant, maxRateConstantChem; // Maximum concentration [mol/m3], diffusion coefficient [m2/s], and homo/hetero max rate constants [1/s] and [??]
public:
    // reaction network structs (vecX is a subset of vecAllX, where vecX contains the items used in the simulation):
    vector<Species*> vecSpecies; // vector holding a s_species struct per species in the system
    vector<Reaction*> vecReactions, vecAllReactions; // vector holding a s_reaction struct per homogeneous reaction in the system
    vector<Redox*> vecRedox, vecAllRedox; // vector holding a s_redox struct per redox couple in the system
    
    System() : maxConcentration( 0.0 ), maxDiffusionConstant( 0.0 ), maxRateConstantChem( 0.0 ) {}

    double getConcMax() const { return maxConcentration; }
    double getDiffMax() const { return maxDiffusionConstant; }
    double getRateMax() const { return maxRateConstantChem; }
    
    vector<double>::size_type addRedox(Redox*);
    vector<double>::size_type addReaction(Reaction*);

    void finalize(double); // parameter: characteristic dimension of system (epsilon)
private:
    void setActiveSystem();
    void setActiveReaction(Reaction*);
    void setActiveRedox(Redox*);
    void setActiveSpecies(Species*);
    template <class sysType> void addToVector(sysType, std::vector<sysType>*);
    bool isActiveSpecies(Species*);

    void updateSpeciesProperties();
    void calcMaxRateConstants(double); // parameter: characteristic dimension of system (epsilon)
    int equilibrateConcentrations();
};

#endif // SYSTEM_H
