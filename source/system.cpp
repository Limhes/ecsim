#include "system.h"
#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

/*===============================================================================================
 * ADD SPECIES, REDOX COUPLES, HOMOGENEOUS REACTIONS TO SYSTEM
 *=============================================================================================*/

Species* System::getSpeciesByName(string name_to_test)
{
    for (auto spec: vecAllSpecies)
    {
        if (spec->getName() == name_to_test) return spec;
    }
    return nullptr;
}

bool System::isSpeciesPresentWithName(string name_to_test)
{
    return (getSpeciesByName(name_to_test) != nullptr);
}

string System::generateUniqueSpeciesName()
{
    // species name will go from A to Z, then into some symbols, then a to z and at one point go wrong (TODO: fix that)
    char c = 'A';
    string s;
    do { s = c++; } while (isSpeciesPresentWithName(s));
    return s;
}

vector<double>::size_type System::addSpecies(Species* spec)
{
    // make sure species has a unique, non-empty name:
    if (spec->getName() == "" || isSpeciesPresentWithName(spec->getName()))
        spec->setName(generateUniqueSpeciesName());

    // add species to system:
    vecAllSpecies.push_back(spec);
    return vecAllSpecies.size();
}

vector<double>::size_type System::addRedox(Redox *redox)
{
    // add redox to system:
    vecAllRedox.push_back(redox);
    return vecAllRedox.size();
}

vector<double>::size_type System::addReaction(Reaction *rxn)
{
    // add reaction to system:
    vecAllReactions.push_back(rxn);
    return vecAllReactions.size();
}

/*===============================================================================================
 * FINALIZE SYSTEM
 *=============================================================================================*/

void System::finalize(double epsilon)
{
    // these functions have to be called in this order:
    setActiveSystem();
    updateSpeciesProperties();
    calcMaxRateConstants(epsilon);
    equilibrateConcentrations();
}

void System::setActiveSystem()
{   
    // clear all vectors, so active redox/reaction/species can be added:
    vecSpecies.clear();
    vecRedox.clear();
    vecReactions.clear();

    // add all redox and contained species:
    for (auto redox: vecAllRedox)
    {
        if (redox->isEnabled())
        {
            // set species to enabled and add to vector that's used by the simulator:
            setActiveSpecies(redox->getSpecRed());
            setActiveSpecies(redox->getSpecOx());
            vecRedox.push_back(redox);
        }
    }

    // add all reaction and contained species:
    for (auto rxn: vecAllReactions)
    {
        if (rxn->isEnabled())
        {
            // set species to enabled and add to vector that's used by the simulator:
            setActiveSpecies(rxn->getSpecLHS1());
            setActiveSpecies(rxn->getSpecLHS2());
            setActiveSpecies(rxn->getSpecRHS1());
            setActiveSpecies(rxn->getSpecRHS2());
            vecReactions.push_back(rxn);
        }
    }
    
    size_t idx = 0;
    for (auto redox: vecRedox) redox->setIndex(idx++);
    idx = 0;
    for (auto rxn: vecReactions) rxn->setIndex(idx++);
    idx = 0;
    for (auto spec: vecSpecies) spec->setIndex(idx++);
}

void System::setActiveSpecies(Species* spec)
{
    // do nothing if "no species" or species has already been activated:
    if (spec == nullptr || isActiveSpecies(spec)) return;
    
    vecSpecies.push_back(spec);
}

bool System::isActiveSpecies(Species* spec_to_test)
{
    for (auto spec: vecSpecies)
    {
        if (spec == spec_to_test) return true;
    }
    return false;
}

void System::updateSpeciesProperties()
{
    // determine maximum concentration and diff coeff:
    maxDiffusionConstant = 0.0;
    maxConcentration = 0.0;
    for (auto spec: vecSpecies)
    {
        maxDiffusionConstant = fmax(maxDiffusionConstant, spec->getDiff());
        maxConcentration = fmax(maxConcentration, spec->getConcInit());
    }
    // normalize concentration and diff coeff, and assign incremental index starting from 0:
    for (auto spec: vecSpecies)
    {
        spec->normalizeDiff(maxDiffusionConstant);
        spec->normalizeConc(maxConcentration);
    }
}

void System::calcMaxRateConstants(double epsilon)
{
    // normalize chemical rates, and find maximum normalized chemical rate constant max_Kchem:
    maxRateConstantChem = 0.0;
    for (auto rxn: vecReactions)
    {
        rxn->normalizeKf(maxDiffusionConstant / (epsilon * epsilon));
        if (rxn->getSpecLHS1() != nullptr && rxn->getSpecLHS2() != nullptr) rxn->normalizeKf(1.0/maxConcentration); // bimolecular rxn
        maxRateConstantChem = fmax(maxRateConstantChem, rxn->getKfNorm());

        rxn->normalizeKb(maxDiffusionConstant / (epsilon * epsilon));
        if (rxn->getSpecRHS1() != nullptr && rxn->getSpecRHS2() != nullptr) rxn->normalizeKb(1.0/maxConcentration); // bimolecular rxn
        maxRateConstantChem = fmax(maxRateConstantChem, rxn->getKbNorm()); // dimensionless
    }
    maxRateConstantChem /= epsilon * epsilon / maxDiffusionConstant; // denormalize for consistency

    // normalize Ke:
    for (auto red: vecRedox)
    {
        red->normalizeKe(maxDiffusionConstant / epsilon); // Ke now dimensionless
    }
}

int System::equilibrateConcentrations()
{
    const double TOO_LONG_TO_WAIT = 1e6;
    const double MIN_CONC_CHANGE = 1e-6; // if (dimensionless) concentration changes by less than this value, we have converged
    const double V_SHRINK = 0.3;
    const double V_GROW = 1.1;

    int underflow, unstable, iter = 0;
    double delta_conc, f, b, minConc1 = 0, minConc2 = 0;
    double dt = 0.1; // since kf, kb, conc are dimensionless!
    
    size_t idx, num_spec = vecSpecies.size();
    vector<double> Del, conc, equilconc;
    
    for (idx = 0; idx < num_spec; idx++)
    {
        Del.push_back(0.0);
        conc.push_back(vecSpecies[idx]->getConcNorm());
        equilconc.push_back(vecSpecies[idx]->getConcNorm());
    }

    do
    {
        for (auto d: Del)
        {
            d = 0.0;
        }

        for (auto rxn: vecReactions)
        {
            iter++;

            f = 1.0;
            b = 1.0;

            for (idx = 0; idx < num_spec; idx++)
            {
                if (rxn->getSpecLHS1() == vecSpecies[idx] || rxn->getSpecLHS2() == vecSpecies[idx]) f *= conc[idx];
                else if (rxn->getSpecRHS1() == vecSpecies[idx] || rxn->getSpecRHS2() == vecSpecies[idx]) b *= conc[idx];
            }

            delta_conc = (f * rxn->getKfNorm() - b * rxn->getKbNorm()) * dt;
            for (idx = 0; idx < num_spec; idx++)
            {
                if (rxn->getSpecLHS1() == vecSpecies[idx] || rxn->getSpecLHS2() == vecSpecies[idx]) Del[idx] -= delta_conc;
                if (rxn->getSpecRHS1() == vecSpecies[idx] || rxn->getSpecRHS2() == vecSpecies[idx]) Del[idx] += delta_conc;
            }
        }

        minConc2 = 1.0;
        underflow = 0;
        unstable = 0;

        for (idx = 0; idx < num_spec; idx++)
        {
            // concentrations equilibrated:
            equilconc[idx] = conc[idx];
            conc[idx] += Del[idx];
            vecSpecies[idx]->setConcNormEquil(equilconc[idx]);

            if (conc[idx] < 0.0)
            {
                underflow++;

                if (conc[idx] < minConc2)
                {
                    minConc2 = conc[idx];
                    minConc1 = equilconc[idx];
                    continue;
                }
            }
            else if (Del[idx] > MIN_CONC_CHANGE)
            {
                unstable++;
            }
        }

        if (underflow > 0)
        {
            dt *= V_SHRINK * minConc1 / (minConc1-minConc2);
            for (idx = 0; idx < num_spec; idx++)
            {
                conc[idx] = equilconc[idx]; // restore original concentrations
            }
        }
        else
        {
            dt *= V_GROW;
        }

    } while (iter < TOO_LONG_TO_WAIT && unstable > 0);

    return iter;
}
