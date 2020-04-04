#include "system.h"
#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

/*===============================================================================================
 * ADD SPECIES, REDOX COUPLES, HOMOGENEOUS REACTIONS TO SYSTEM
 *=============================================================================================*/

bool System::isSpeciesPresentWithName(string _name)
{
    for (auto spec: vecAllSpecies)
    {
        if (spec->name == _name) return true;
    }
    return false;
}

Species* System::getSpeciesByName(string _name)
{
    for (auto spec: vecAllSpecies)
    {
        if (spec->name == _name) return spec;
    }
    return nullptr;
}

string System::generateUniqueSpeciesName()
{
    // species name will go from A to Z, then into some symbols, then a to z and at one point go wrong...
    char c = 'A';
    string s;
    do { s = c++; } while (isSpeciesPresentWithName(s));
    return s;
}

vector<double>::size_type System::addSpecies(Species *spec)
{
    if (spec->name == "") spec->name = generateUniqueSpeciesName();

    // add species to system:
    vecAllSpecies.emplace_back(spec);

    return vecAllSpecies.size();
}

vector<double>::size_type System::addRedox(Redox *redox)
{
    // add redox to system:
    vecAllRedox.emplace_back(redox);

    return vecAllRedox.size();
}

vector<double>::size_type System::addReaction(Reaction *rxn)
{
    // add reaction to system:
    vecAllReactions.emplace_back(rxn);

    return vecAllReactions.size();
}

/*===============================================================================================
 * FINALIZE SYSTEM
 *=============================================================================================*/

void System::finalize(double epsilon)
{
    setActiveSystem();
    updateSpeciesProperties();
    calcMaxRateConstants(epsilon);
    equilibrateConcentrations();
}

void System::setActiveSystem()
{
    // first, set all species to disabled:
    for (auto spec: vecAllSpecies)
    {
        spec->inRedox = false;
        spec->inReaction = false;
    }

    // add all redox:
    vecRedox.clear();
    size_t redoxidx = 0;
    for (auto redox: vecAllRedox)
    {
        if (redox->enabled)
        {
            // assign index, set species to enabled and add to vector that's used by the simulator:
            redox->index = redoxidx++;
            if (redox->specReduced != nullptr)  redox->specReduced->inRedox = true;
            if (redox->specOxidized != nullptr) redox->specOxidized->inRedox = true;
            vecRedox.emplace_back(redox);
        }
    }

    // add all reaction:
    vecReactions.clear();
    size_t rxnidx = 0;
    for (auto rxn: vecAllReactions)
    {
        if (rxn->enabled)
        {
            // assign index, set species to enabled and add to vector that's used by the simulator:
            rxn->index = rxnidx++;
            if (rxn->specLHS1 != nullptr) rxn->specLHS1->inReaction = true;
            if (rxn->specLHS2 != nullptr) rxn->specLHS2->inReaction = true;
            if (rxn->specRHS1 != nullptr) rxn->specRHS1->inReaction = true;
            if (rxn->specRHS2 != nullptr) rxn->specRHS2->inReaction = true;
            vecReactions.emplace_back(rxn);
        }
    }

    // add all species:
    vecSpecies.clear();
    size_t specidx = 0;
    for (auto spec: vecAllSpecies)
    {
        if (spec->inRedox || spec->inReaction)
        {
            // assign index and add to vector that's used by the simulator:
            spec->index = specidx++;
            vecSpecies.emplace_back(spec);
        }
    }
}

void System::updateSpeciesProperties()
{
    // determine maximum concentration and diff coeff:
    maxDiffusionConstant = 0.0;
    maxConcentration = 0.0;
    for (auto spec: vecSpecies)
    {
        maxDiffusionConstant = fmax(maxDiffusionConstant, spec->diffusionConstant);
        maxConcentration = fmax(maxConcentration, spec->initialConcentration);
    }
    // normalize concentration and diff coeff, and assign incremental index starting from 0:
    for (auto spec: vecSpecies)
    {
        spec->normalizedDiffusionConstant = spec->diffusionConstant / maxDiffusionConstant;
        spec->normalizedConcentration = spec->initialConcentration / maxConcentration;
    }
}

void System::calcMaxRateConstants(double epsilon)
{
    // normalize chemical rates, and find maximum normalized chemical rate constant max_Kchem:
    maxRateConstantChem = 0.0;
    for (auto rxn: vecReactions)
    {
        rxn->rateConstantForwardNormalized = rxn->rateConstantForward * epsilon * epsilon / maxDiffusionConstant;
        if (rxn->specLHS1 != nullptr && rxn->specLHS2 != nullptr) rxn->rateConstantForwardNormalized *= maxConcentration; // bimolecular rxn
        maxRateConstantChem = fmax(maxRateConstantChem, rxn->rateConstantForwardNormalized);

        rxn->rateConstantBackwardNormalized = rxn->rateConstantBackward * epsilon * epsilon / maxDiffusionConstant;
        if (rxn->specRHS1 != nullptr && rxn->specRHS2 != nullptr) rxn->rateConstantBackwardNormalized *= maxConcentration; // bimolecular rxn
        maxRateConstantChem = fmax(maxRateConstantChem, rxn->rateConstantBackwardNormalized); // dimensionless
    }
    maxRateConstantChem /= epsilon * epsilon / maxDiffusionConstant; // denormalize for consistency

    // normalize Ke:
    for (auto red: vecRedox)
    {
        red->rateConstantHeteroNormalized = red->rateConstantHetero * epsilon / maxDiffusionConstant; // Ke now dimensionless
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

    for (auto spec: vecSpecies)
    {
        spec->Del = 0.0;
        spec->conc = spec->normalizedConcentration;
        spec->normalizedAndEquilibratedConcentration = spec->normalizedConcentration;
    }

    do
    {
        for (auto spec: vecSpecies)
        {
            spec->Del = 0.0;
        }

        for (auto rxn: vecReactions)
        {
            iter++;

            f = 1.0;
            b = 1.0;

            for (auto spec: vecSpecies)
            {
                if (rxn->specLHS1 == spec || rxn->specLHS2 == spec) f *= spec->conc;
                else if (rxn->specRHS1 == spec || rxn->specRHS2 == spec) b *= spec->conc;
            }

            delta_conc = (f * rxn->rateConstantForwardNormalized - b * rxn->rateConstantBackwardNormalized) * dt;
            for (auto spec: vecSpecies) {
                if (rxn->specLHS1 == spec || rxn->specLHS2 == spec) spec->Del -= delta_conc;
                if (rxn->specRHS1 == spec || rxn->specRHS2 == spec) spec->Del += delta_conc;
            }
        }

        minConc2 = 1.0;
        underflow = 0;
        unstable = 0;

        for (auto spec: vecSpecies)
        {
            spec->normalizedAndEquilibratedConcentration = spec->conc;
            spec->conc += spec->Del;

            if (spec->conc < 0.0)
            {
                underflow++;

                if (spec->conc < minConc2)
                {
                    minConc2 = spec->conc;
                    minConc1 = spec->normalizedAndEquilibratedConcentration;
                    continue;
                }
            }
            else if (spec->Del > MIN_CONC_CHANGE)
            {
                unstable++;
            }
        }

        if (underflow > 0)
        {
            dt *= V_SHRINK * minConc1 / (minConc1-minConc2);
            for (auto spec: vecSpecies)
            {
                spec->conc = spec->normalizedAndEquilibratedConcentration; // restore original concentrations
            }
        }
        else
        {
            dt *= V_GROW;
        }

    } while (iter < TOO_LONG_TO_WAIT && unstable > 0);

    return iter;
}


/*===============================================================================================
 * DESTRUCT ALL
 *=============================================================================================*/

System::~System()
{
    vecAllRedox.clear();
    vecAllReactions.clear();
    vecAllSpecies.clear();
}
