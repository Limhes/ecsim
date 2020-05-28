#include "system.h"
#include <cmath>
#include <iostream>
//#include <algorithm>
#include <array>

/*===============================================================================================
 * ADD REDOX STEPS, HOMOGENEOUS REACTIONS TO SYSTEM
 *=============================================================================================*/

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
    setActiveSystem(); // fills vecSpecies/vecRedox/vecReactions with the active moieties
    updateSpeciesProperties(); // fills maxConcentration/maxDiffusionConstant and normalizes species conc/diff
    calcMaxRateConstants(epsilon); // fills maxRateConstantChem and normalized reactions kf/kb
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
        setActiveRedox(redox);
    }

    // add all reaction and contained species:
    for (auto rxn: vecAllReactions)
    {
        setActiveReaction(rxn);
    }
}

void System::setActiveReaction(Reaction* rxn)
{
    if (rxn->isEnabled())
    {
        // set species to enabled and add to vector that's used by the simulator:
        setActiveSpecies(rxn->getSpecLHS1());
        setActiveSpecies(rxn->getSpecLHS2());
        setActiveSpecies(rxn->getSpecRHS1());
        setActiveSpecies(rxn->getSpecRHS2());
        
        addToVector(rxn, &vecReactions);
    }
}

void System::setActiveRedox(Redox* redox)
{
    if (redox->isEnabled())
    {
        // set species to enabled and add to vector that's used by the simulator:
        setActiveSpecies(redox->getSpecRed());
        setActiveSpecies(redox->getSpecOx());
        
        addToVector(redox, &vecRedox);
    }
}

void System::setActiveSpecies(Species* spec)
{
    // do nothing if "no species" or species has already been activated:
    if (spec != nullptr && !isActiveSpecies(spec))
    {
        addToVector(spec, &vecSpecies);
    }
}

template <class sysType>
void System::addToVector(sysType obj, std::vector<sysType>* vec)
{
    obj->setIndex(vec->size());
    vec->push_back(obj);
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
        maxDiffusionConstant = std::fmax(maxDiffusionConstant, spec->getDiff());
        maxConcentration = std::fmax(maxConcentration, spec->getConcInit());
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
    double normfactor;
    // normalize chemical rates, and find maximum normalized chemical rate constant max_Kchem:
    maxRateConstantChem = 0.0;
    for (auto rxn: vecReactions)
    {
        normfactor = maxDiffusionConstant / (epsilon * epsilon);
        if (rxn->getSpecLHS1() != nullptr && rxn->getSpecLHS2() != nullptr) normfactor /= maxConcentration; // bimolecular rxn
        rxn->normalizeKf(normfactor);
        maxRateConstantChem = std::fmax(maxRateConstantChem, rxn->getKfNorm());

        normfactor = maxDiffusionConstant / (epsilon * epsilon);
        if (rxn->getSpecRHS1() != nullptr && rxn->getSpecRHS2() != nullptr) normfactor /= maxConcentration; // bimolecular rxn
        rxn->normalizeKb(normfactor);
        maxRateConstantChem = std::fmax(maxRateConstantChem, rxn->getKbNorm()); // dimensionless
    }
    maxRateConstantChem /= epsilon * epsilon / maxDiffusionConstant; // denormalize for consistency

    // normalize Ke:
    for (auto red: vecRedox)
    {
        red->normalizeKe(maxDiffusionConstant / epsilon); // dimensionless
    }
}

int System::equilibrateConcentrations()
{
    /*
     * since we do not know what type of system we have to equilibrate, we have to simulate it using forward and backward reactions
     * this is an implementation of Bray, ...
     */
    const double TOO_LONG_TO_WAIT = 1e6;
    const double MIN_CONC_CHANGE = 1e-6; // if (dimensionless) concentration changes by less than this value, we have converged
    const double V_SHRINK = 0.3;
    const double V_GROW = 1.1;

    int underflow, unstable, iter = 0;
    double delta_conc, f, b, minConc1 = 0.0, minConc2 = 0.0;
    double dt = 0.1; // kf, kb, conc are dimensionless!
    
    size_t idx, num_spec = vecSpecies.size();
    std::vector<double> Del(num_spec), conc(num_spec), equilconc(num_spec);
    
    for (idx = 0; idx < num_spec; idx++)
    {
        conc[idx] = vecSpecies[idx]->getConcNorm();
        equilconc[idx] = vecSpecies[idx]->getConcNorm();
    }

    do
    {
        iter++;
        
        for (idx = 0; idx < num_spec; idx++)
        {
            Del[idx] = 0.0;
        }

        for (auto rxn: vecReactions)
        {
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
