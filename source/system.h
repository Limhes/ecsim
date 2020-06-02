#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include <iostream>
#include <cassert>

const std::string NO_SPECIES_NAME = "(no species)";

class Species {
private:
    double initialConcentration, normalizedConcentration, normalizedAndEquilibratedConcentration;
    double diffusionConstant, normalizedDiffusionConstant;
    size_t index;
    std::string name;
public:
    double getConcInit() const { return initialConcentration; }
    void setConcInit(double c) { initialConcentration = c; }
    double getConcEquil() const { assert(normalizedConcentration > 0.0); return normalizedAndEquilibratedConcentration*initialConcentration/normalizedConcentration; }
    double getDiff() const { return diffusionConstant; }
    void setDiff(double d) { diffusionConstant = d; }
    
    double getConcNorm() const { return normalizedConcentration; }
    double getConcNormEquil() const { return normalizedAndEquilibratedConcentration; }
    void setConcNormEquil(double c) { normalizedAndEquilibratedConcentration = c; }
    double getDiffNorm() const { return normalizedDiffusionConstant; }
    void normalizeConc(double conc_max) { assert(conc_max > 0.0); normalizedConcentration = initialConcentration/conc_max; }
    void normalizeDiff(double diff_max) { assert(diff_max > 0.0); normalizedDiffusionConstant = diffusionConstant/diff_max; }
    
    void setName(std::string new_name) { name = new_name; }
    std::string getName() const { return name; }
    void setIndex(std::size_t new_index) { index = new_index; }
    std::size_t getIndex() const { return index; }

    Species(std::string n, double c, double d)
      : initialConcentration( c ), normalizedConcentration( 0.0 ), normalizedAndEquilibratedConcentration( 0.0 ),
        diffusionConstant( d ), normalizedDiffusionConstant( 0.0 ), index( 0 ), name( n )
    {}
};

class Reaction
{ // specLHS1 + specLHS2 <==> specRHS1 + specRHS2
    Species *specLHS1, *specLHS2, *specRHS1, *specRHS2;
    double rateConstantForward, rateConstantForwardNormalized; // holds the normalized forward rate constant: Kf * t (unimolecular) or Kf * t * max_conc (bimolecular)
    double rateConstantBackward, rateConstantBackwardNormalized; // holds the normalized backward rate constant: Kb * t (unimolecular) or Kb * t * max_conc (bimolecular)
    std::size_t index;
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
    void normalizeKf(double norm) { assert(norm > 0.0); rateConstantForwardNormalized = rateConstantForward / norm; }
    void normalizeKb(double norm) { assert(norm > 0.0); rateConstantBackwardNormalized = rateConstantBackward / norm; }
    
    void setIndex(std::size_t new_index) { index = new_index; }
    std::size_t getIndex() const { return index; }
    
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
    unsigned int numberElectrons; // holds the number of electrons
    double standardPotential; // holds the potential [V]
    double rateConstantHetero, rateConstantHeteroNormalized; // holds the normalized heterogeneous rate constant: Ke * delta_t / delta_x
    double alpha; // holds the alpha coefficient
    std::size_t index;
    bool enabled;
public:
    Species* getSpecOx() const { return specOxidized; }
    Species* getSpecRed() const { return specReduced; }
    
    double getKe() const { return rateConstantHetero; }
    void setKe(double k) { rateConstantHetero = k; }
    double getKeNorm() const { return rateConstantHeteroNormalized; }
    void normalizeKe(double norm) { assert(norm > 0.0); rateConstantHeteroNormalized = rateConstantHetero / norm; }
    double getE0() const { return standardPotential; }
    void setE0(double E) { standardPotential = E; }
    unsigned int getNe() const { return numberElectrons; }
    void setNe(unsigned int n) { numberElectrons = n; }
    double getAlpha() const { return alpha; }
    void setAlpha(double a) { alpha = a; }
    
    void setIndex(std::size_t new_index) { index = new_index; }
    std::size_t getIndex() const { return index; }
    
    bool isEnabled() const { return enabled; }
    void enable() { enabled = true; }
    void disable() { enabled = false; }

    Redox(Species* ox, Species* red, unsigned int n, double E, double Ke, double a)
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
    std::vector<Species*> vecSpecies; // vector holding a s_species struct per species in the system
    std::vector<Reaction*> vecReactions, vecAllReactions; // vector holding a s_reaction struct per homogeneous reaction in the system
    std::vector<Redox*> vecRedox, vecAllRedox; // vector holding a s_redox struct per redox couple in the system
    
    System() : maxConcentration( 0.0 ), maxDiffusionConstant( 0.0 ), maxRateConstantChem( 0.0 ) {}

    double getConcMax() const { return maxConcentration; }
    double getDiffMax() const { return maxDiffusionConstant; }
    double getRateMax() const { return maxRateConstantChem; }
    
    std::vector<double>::size_type addRedox(Redox*);
    std::vector<double>::size_type addReaction(Reaction*);

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
    unsigned int equilibrateConcentrations();
};

#endif // SYSTEM_H
