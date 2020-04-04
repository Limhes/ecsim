#include <iostream>
#include "simulation.h"

using namespace std;

uint64_t getPosixClockTime()
{
    struct timespec ts = {};
    uint64_t elapsed_time_ns;

    if (clock_gettime (CLOCK_MONOTONIC, &ts) == 0)
    {
        elapsed_time_ns = static_cast<uint64_t>(ts.tv_sec * 1e9 + ts.tv_nsec);
    }
    else
    {
        elapsed_time_ns = 0;
    }

    return elapsed_time_ns;
}

/*===============================================================================================
 * HIGHEST LEVEL CLASS METHODS
 *=============================================================================================*/

size_t Simulation::run(vector<double> &_current, vector<double> &_potential)
{
    uint64_t startTime, endTime; // timing variables [ns]

    // prepare simulation:
    prepareSimulation();

    // run simulation
    startTime = getPosixClockTime();
    size_t numDataPoints = runSimulation(_current, _potential);
    endTime = getPosixClockTime();
    output_stream << "Est. simulation time (" << numDataPoints << " steps &amp; " << numGridPoints << " grid points) = " << static_cast<double>(endTime - startTime)/1.0e9 << " s<br />" << endl;

    return numDataPoints;
}

size_t Simulation::runSimulation(vector<double> &_current, vector<double> &_potential)
{
    double segmentStartPotential; // start potential of scan segment
    bool recordCurrent; // record current in segment?

    matrixPatternAnalyzed = false;

    // initialize electroanalytical parameters:
    ipc = 0.0;
    ipa = 0.0;
    Epc = 0.0;
    Epa = 0.0;

    // conditioning step:
    delaySegment(exper->conditioningPotential, exper->conditioningTime);
    // equilibration step:
    delaySegment(exper->initialPotential, exper->equilibrationTime);
    // scan cycles:
    for (int cycle = 0; cycle < exper->numCycles; cycle++) {
        recordCurrent = (cycle == exper->numCycles - 1); // warning: hard-coded option to record current only of last cycle

        segmentStartPotential = exper->initialPotential;
        for (auto vp: exper->vertexPotentials) {
            // scan up to vertex:
            scanSegment(segmentStartPotential, vp, recordCurrent, _current, _potential);
            segmentStartPotential = vp;
            // vertex delay:
            delaySegment(vp, exper->vertexDelay);
        }
        // scan from last vertex to final potential:
        scanSegment(segmentStartPotential, exper->finalPotential, recordCurrent, _current, _potential);
    }

    return _current.size();
}

void Simulation::delaySegment(double potential, double time)
{
    if (time > MIN_STEP_TIME)
    {
        double theta = f*potential;

        int num_points = static_cast<int>(time/(totalTime*deltaTheta/totalTheta));
        for (int d = 0; d < num_points; d++)
        {
            // do magic:
            solveSystem(theta);
        }
    }
}

void Simulation::scanSegment(double potential_start, double potential_stop, bool record_current, vector<double> &_current, vector<double> &_potential)
{
    double curr;
    double theta_start = f*potential_start;
    double theta_stop = f*potential_stop;

    double sign = 0.0;
    if (theta_start < theta_stop) sign = 1.0;
    else if (theta_start > theta_stop) sign = -1.0;

    int num_points = static_cast<int>(abs(theta_start-theta_stop)/deltaTheta);
    double theta = theta_start;
    for (int d = 0; d < num_points; d++)
    {
        // do magic:
        solveSystem(theta);

        if (record_current) {
            // get current:
            curr = calcCurrentFromFlux();

            // electroanalytical metrics:
            if (curr > ipa) { ipa = curr; Epa = theta/f; }
            if (curr < ipc) { ipc = curr; Epc = theta/f; }

            // add data points to current & potential vectors:
            _potential.emplace_back(theta/f);
            _current.emplace_back(curr);
        }

        // increase potential:
        theta += sign*deltaTheta;
    }
}

/*===============================================================================================
 * THIS IS THE SIMULATION
 *=============================================================================================*/

void Simulation::solveSystem(double theta)
{
    storeConcentrations();
    updateIndependentTerms();
    updateRedoxInMatrix(theta); // update potentials in matrix
    updateKineticsInMatrix(true);
    if (!matrixPatternAnalyzed) { sparseMatrixSolver.analyzePattern(matrixA); matrixPatternAnalyzed = true; }
    invertMatrix(); // Ax=b, solve for x
    updateKineticsInMatrix(false);
}

void Simulation::storeConcentrations()
{
    for (auto spec: sys->vecSpecies)
    {
        for (size_t x = 0; x < numGridPoints; x++)
        {
            gridConcentrationPrevious[spec->index*numGridPoints+x] = gridConcentration[spec->index*numGridPoints+x];
        }
    }
}

void Simulation::updateKineticsInMatrix(bool add)
{
    double rate;
    size_t spec1, spec2, spec3;
    Eigen::Index idx1, idx2, idx3;

    for (auto kinTerm: secondOrderKinetics)
    {
        spec1 = get<0>(kinTerm);
        spec2 = get<1>(kinTerm);
        spec3 = get<2>(kinTerm);

        for (size_t x = 1; x < numGridPoints; x++)
        {
            rate = get<3>(kinTerm) * deltaX*deltaX*paramGamma2i[x];
            idx1 = static_cast<Eigen::Index>(spec1*numGridPoints+x);
            idx2 = static_cast<Eigen::Index>(spec2*numGridPoints+x);
            idx3 = static_cast<Eigen::Index>(spec3*numGridPoints+x);

            if (add)
            {
                matrixA.coeffRef( idx3, idx1 ) += rate * gridConcentration[spec2*numGridPoints+x];
                matrixA.coeffRef( idx3, idx2 ) += rate * gridConcentration[spec1*numGridPoints+x];
                independentTerms[spec3*numGridPoints+x] += rate * gridConcentration[spec1*numGridPoints+x] * gridConcentration[spec2*numGridPoints+x];
            }
            else
            {
                matrixA.coeffRef( idx3, idx1 ) -= rate * gridConcentrationPrevious[spec2*numGridPoints+x];
                matrixA.coeffRef( idx3, idx2 ) -= rate * gridConcentrationPrevious[spec1*numGridPoints+x];
            }
        }
    }
}

void Simulation::updateIndependentTerms()
{
    size_t EndNormal = numGridPoints - numDerivPoints + 2;

    //for(int s = 0; s < sys->num_species; s++)
    for (auto spec: sys->vecSpecies)
    {
        //Zero at matrix rows storing surface conditions
        independentTerms[spec->index*numGridPoints] = 0;

        for (size_t x = 1; x < numGridPoints; x++)
        {
            /* independentTerms is the RHS ('b') of the matrix equation A*x=b, where 'x' is the concentration at T+dT (the 'next' concentration)
             * the system of equations is ('a' are coefficients, 'x(i,s)' the next concentrations at grid point i for species s):
             * a(-1)*x(i-1,s) + a(0)*x(i,s) + a(+1)*x(i+1,s) + a(+2)*x(i+2,s) + (etc.) = b(i,s)
             */
            independentTerms[x+spec->index*numGridPoints] = -gridConcentration[x+numGridPoints*spec->index]*paramGamma2i[x]/paramLambda;
            // however, when i>numGridPoints, then x(i+i) reduces to the (never changing) bulk concentration and we add them into 'b', because they are now known:
            if (x >= EndNormal)
            {
                for (int jj = 0; jj <= static_cast<int>(x - EndNormal); jj++)
                {
                    independentTerms[x+spec->index*numGridPoints] -= coeffMatrixN2(x, static_cast<int>(numDerivPoints)-jj-2, spec->normalizedDiffusionConstant) * spec->normalizedAndEquilibratedConcentration;
                }
            }
        }
    }
}

void Simulation::updateRedoxInMatrix(double theta)
{
    Eigen::Index oxidx, redidx, specidx;
    double p, Kred, Kox;

    // reset diagonal terms:
    fill(matrixB0DiagonalTerms.begin(), matrixB0DiagonalTerms.end(), 0.0);

    for (auto redox: sys->vecRedox)
    {
        p = static_cast<double>(redox->numberElectrons) * (theta - f * redox->standardPotential); // normalized potential

        Kred = redox->rateConstantHeteroNormalized * exp(-redox->alpha * p); // B-V kinetics
        Kox = redox->rateConstantHeteroNormalized * exp((1.0 - redox->alpha) * p); // B-V kinetics

        oxidx = static_cast<Eigen::Index>(redox->specOxidized->index*numGridPoints);
        redidx = static_cast<Eigen::Index>(redox->specReduced->index*numGridPoints);

        // add terms per redox reaction:
        matrixA.coeffRef( oxidx,  redidx ) = - deltaX * Kox / redox->specOxidized->normalizedDiffusionConstant;  // B0
        matrixA.coeffRef( redidx, oxidx )  = - deltaX * Kred / redox->specReduced->normalizedDiffusionConstant; // B0

        // add to diagonal terms:
        matrixB0DiagonalTerms[redox->specOxidized->index] += deltaX * Kred / redox->specOxidized->normalizedDiffusionConstant;
        matrixB0DiagonalTerms[redox->specReduced->index] += deltaX * Kox / redox->specReduced->normalizedDiffusionConstant;
    }

    // add diagonal terms:
    for (auto spec: sys->vecSpecies)
    {
        if (speciesInRedox[spec->index]) // active redox species
        {
            specidx = static_cast<int>(spec->index*numGridPoints);
            matrixA.coeffRef( specidx, specidx )   = 1.0 + matrixB0DiagonalTerms[spec->index]; // B0
            matrixA.coeffRef( specidx, specidx+1 ) = -1.0; // B1
        }
    }
}

void Simulation::invertMatrix()
{
    sparseMatrixSolver.factorize(matrixA);

    for (size_t kk = 0; kk < numOneRow; kk++) vecb[static_cast<Eigen::Index>(kk)] = independentTerms[kk];
    vecx = sparseMatrixSolver.solve(vecb);
    for (size_t kk = 0; kk < numOneRow; kk++) gridConcentration[kk] = vecx[static_cast<Eigen::Index>(kk)];
}

double Simulation::calcCurrentFromFlux()
{
    double totalflux, speciesflux, currentFromFlux;
    currentFromFlux = el->epsilon * CONST_F * sys->maxDiffusionConstant * sys->maxConcentration;

    // determine flux per species:
    for (auto spec: sys->vecSpecies)
    {
        speciesflux = 0;
        for(size_t kk = 0; kk < numCurrentPoints; kk++)
        {
            speciesflux += coeffBeta0[kk] * gridConcentration[spec->index*numGridPoints+kk];
        }
        currentContributionSpeciesFlux[static_cast<Eigen::Index>(spec->index)] = speciesflux * spec->normalizedDiffusionConstant;
    }

    /* SOLVE FOR THE FLUX PER REDOX REACTION BY LEAST SQUARES
     * - The method is described in Compton (pages 117 through 119) but is not general
     * - I have generalized the method as follows, and it works for all mechanisms I have tried:
     *   - The nett flux per redox reaction flux_redox = k_red*[ox] - k_ox*[red]
     *   - The total flux at the electrode (from which the current can be calculated) can be expressed as:
     *     total_flux = - SUM(flux_redox, for each redox reaction)
     *   - The flux of each species at the electrode is:
     *     +redox_flux if the species is ox in that redox reaction
     *     -redox_flux if the species is red in that redox reaction
     *     0 if it doesn't participate in that redox reaction
     *     --> I have put all these "current contributions" in currentContributionMatrix(species, redox) where each value is -1.0/0.0/1.0
     *   - The flux of each species at the electrode is ALSO:
     *     - Calculated by the "current function", using beta0 values (see for loop above)
     *     --> These species fluxes are stored in currentContributionSpeciesFlux(species)
     *   - The nett flux per redox reaction is then determined by solving the (overdetermined!) matrix equation:
     *     F*fr=fs where F = currentContributionMatrix, fr = flux per redox reaction, fs = species flux
     *   - The total flux is then determined by summing the elements in fr
     *   :)
     */
    // determine nett flux per redox reaction
    currentContributionRedoxFlux = (currentContributionMatrix.transpose() * currentContributionMatrix).ldlt().solve(
                currentContributionMatrix.transpose() * currentContributionSpeciesFlux );

    // sum to get total flux:
    totalflux = 0.0;
    for (auto redox: sys->vecRedox)
    {
        // this doesn't work well for n > 1 --> why not?
        totalflux -= currentContributionRedoxFlux[static_cast<Eigen::Index>(redox->index)]
                     * static_cast<double>(redox->numberElectrons);
    }
    totalflux /= deltaX; // totalflux := dC/dX

    // calculate current from flux:
    return totalflux * currentFromFlux;
}

/*===============================================================================================
 * BEFORE SIMULATION
 *=============================================================================================*/

void Simulation::prepareSimulation()
{
    // initialize all simulation/grid parameters:
    initParametersEtc();
    // initialize all vectors and matrices:
    initVectorsEtc();

    // set boundary conditions in matrix:
    addRedoxToMatrix();
    // fill matrix with backwards implicit (Laasonen) coefficients:
    addBICoeffsToMatrix();
    // add (1st and 2nd order) kinetics to matrix:
    addKineticsToMatrix();
    // create matrix:
    createMatrix();
}

void Simulation::initParametersEtc()
{
    // some parameters
    sys->finalize(el->epsilon);
    f = CONST_F / (CONST_R * env->temperature);
    totalTime = exper->totalTime();
    totalTheta = f * totalTime * exper->scanRate;
    paramR0 = el->epsilon / sqrt(sys->maxDiffusionConstant*totalTime);

    /* system dimensioning according to Compton (Understanding Voltammetry)
     * deltaTheta is a fixed value, taken from Compton
     * to increase accuracy, we need to decrease deltaX, which doesn't need a lot of extra time to solve!!!
     * the factor F in deltaX = pow(10.0, -F-0.5*log10(sigma)); should e.g. go from 2.2 to 4.2
     */
    double F,logRate = log10(sys->maxRateConstantChem);
    if (logRate < minLogRate) { F = minF; }
    else if (logRate <= maxLogRate) { F = minF + (logRate - minLogRate) * (maxF - minF) / (maxLogRate - minLogRate); }
    else { F = maxF; }

    sigma = el->epsilon * el->epsilon / sys->maxDiffusionConstant * f * exper->scanRate; // dimensionless scan rate
    maxT = totalTheta / sigma;
    deltaT = deltaTheta / sigma;
    maxX = 6.0*sqrt(maxT);
    //deltaX = pow(10.0, -deltaXSizingFactor-0.5*log10(sigma));
    deltaX = pow(10.0, -F) / sqrt(sigma); // Compton (page 64) leads to: dX = 10^(-2.2-0.5*sigma)
    paramLambda = deltaT / (deltaX * deltaX);

    // print scaling parameters:
    //output_stream << "Diff coeff [m2/s]: max = " << sys->maxDiffusionConstant << endl;
    //output_stream << "Conc [mol/m3]: max = " << sys->maxConcentration << endl;
    //output_stream << "Electrode [m]: epsilon = " << el->epsilon << ", area = " << el->electrodeArea << endl;
    output_stream << "Initial concentrations after equilibration: (";
    for (auto spec: sys->vecSpecies)
        output_stream << "[" << spec->name << "] = " << spec->normalizedAndEquilibratedConcentration << ((spec!=sys->vecSpecies.back())?", ":"");
    output_stream << ")<br />" << endl;
    output_stream << "V<sub>step</sub> [V] = " << deltaTheta/f << ", t<sub>max</sub> [s] = " << totalTime << ", " << endl;
    output_stream << "r<sub>max</sub> [1/s] = " << sys->maxRateConstantChem << ", F = " << F << "<br />" << endl;
    output_stream << "&sigma; = " << sigma << ", &theta;<sub>max</sub> = " << totalTheta << ", &Delta;&theta; = " << deltaTheta << ", " << endl;
    output_stream << "X<sub>max</sub> = " << maxX << ", &Delta;X = " << deltaX << ", T<sub>max</sub> = " << maxT << ", &Delta;T = " << deltaT << "<br />" << endl;

    // create expanding grid and set dimensions:
    numGridPoints = 1;
    do { numGridPoints++; } while (deltaX < maxX * (paramGamma - 1.0) / ( pow(paramGamma, numGridPoints-1) - 1.0 ));
    numCurrentPoints = min(numCurrentPoints, numGridPoints); // truncate number of current points if it is larger than number of grid points
    numOneRow = numGridPoints * sys->vecSpecies.size(); // == MyData.OneRow

    // output parametrisation of grid:
    //output_stream << "Number of grid points = " << numGridPoints << "<br />" << endl;
    //output_stream << "Number of current points = " << numCurrentPoints << endl;
    //output_stream << "Number of deriv points = " << numDerivPoints << endl;
}

void Simulation::initVectorsEtc()
{
    size_t numNonZeroesInRow = numDerivPoints*2; // approximate value (in excess) of Non Zeros per row
    size_t numNonZeroesTotal = numNonZeroesInRow*numOneRow;

    // initialize distance in grid, expansion factor gamma and concentration vectors:
    gridCoordinate.resize(numGridPoints);
    paramGammai.resize(numGridPoints);
    paramGamma2i.resize(numGridPoints);
    gridConcentration.resize(numOneRow);
    gridConcentrationPrevious.resize(numOneRow);
    for (size_t x = 0; x < numGridPoints; x++)
    {
        paramGammai[x] = pow(paramGamma, x);
        paramGamma2i[x] = paramGammai[x]*paramGammai[x];
        gridCoordinate[x] = (x == 0) ? 0.0 : gridCoordinate[x-1] + deltaX*paramGammai[x-1];

        for (auto spec: sys->vecSpecies)
        {
            gridConcentration[spec->index*numGridPoints+x] = spec->normalizedAndEquilibratedConcentration;
            gridConcentrationPrevious[spec->index*numGridPoints+x] = 0.0;
        }
    }

    // initialize independent terms and b & x vectors:
    independentTerms.resize(numOneRow);
    vecb.resize(static_cast<Eigen::Index>(numOneRow));
    vecx.resize(static_cast<Eigen::Index>(numOneRow));
    for (size_t n = 0; n < numOneRow; n++)
    {
        independentTerms[n] = 0.0;
        vecb[static_cast<Eigen::Index>(n)] = 0.0;
        vecx[static_cast<Eigen::Index>(n)] = 0.0;
    }

    // obtain flux condition coefficients:
    coeffBeta0.resize(numCurrentPoints);
    for(size_t i = 0; i < numCurrentPoints; i++)
    {
        coeffBeta0[i] = Beta_N_1(static_cast<int>(numCurrentPoints), static_cast<int>(i), paramGamma);
    }
    // obtain diffusion coefficients:
    coeffAlpha.resize(numDerivPoints);
    coeffBeta.resize(numDerivPoints);
    for(size_t d = 0; d < numDerivPoints; d++)
    {
        coeffAlpha[d] = Alpha_N_2(static_cast<int>(numDerivPoints), static_cast<int>(d)-1, paramGamma);
        coeffBeta[d] = Beta_N_2(static_cast<int>(numDerivPoints), static_cast<int>(d)-1, paramGamma);
    }

    // initiliaze redox vectors & matrix:
    speciesInRedox.clear();
    speciesInRedox.resize(sys->vecSpecies.size(), false);
    matrixB0DiagonalTerms.clear();
    matrixB0DiagonalTerms.resize(sys->vecSpecies.size(), 0.0);
    currentContributionMatrix.resize(static_cast<Eigen::Index>(sys->vecSpecies.size()), static_cast<Eigen::Index>(sys->vecRedox.size()));
    currentContributionMatrix.setZero();
    currentContributionSpeciesFlux.resize(static_cast<Eigen::Index>(sys->vecSpecies.size()));
    currentContributionRedoxFlux.resize(static_cast<Eigen::Index>(sys->vecRedox.size()));

    // initialize Triplet container
    tripletContainerMatrixA.clear();
    tripletContainerMatrixA.reserve(numNonZeroesTotal);

    // clear and resize matrix:
    matrixA.setZero(); // removes all non-zeros
    matrixA.resize(0,0);
    matrixA.data().squeeze(); // removes all items from matrix
    matrixA.resize(static_cast<Eigen::Index>(numOneRow), static_cast<Eigen::Index>(numOneRow));
    matrixA.reserve(static_cast<Eigen::Index>(numNonZeroesTotal));
}

void Simulation::addKineticsToMatrix()
{
    Species *s;
    for (auto rxn: sys->vecReactions) {
        // FORWARD REACTIONS
        if (rxn->specLHS1 != nullptr && rxn->specLHS2 != nullptr)
        { // second-order forward reaction
            addSecondOrderKineticTerm(rxn->specLHS1->index, rxn->specLHS2->index, rxn->specLHS1->index, -(rxn->rateConstantForwardNormalized));
            addSecondOrderKineticTerm(rxn->specLHS1->index, rxn->specLHS2->index, rxn->specLHS2->index, -(rxn->rateConstantForwardNormalized));
            if (rxn->specRHS1 != nullptr) addSecondOrderKineticTerm(rxn->specLHS1->index, rxn->specLHS2->index, rxn->specRHS1->index, rxn->rateConstantForwardNormalized);
            if (rxn->specRHS2 != nullptr) addSecondOrderKineticTerm(rxn->specLHS1->index, rxn->specLHS2->index, rxn->specRHS2->index, rxn->rateConstantForwardNormalized);
        }
        else
        { // first-order forward reaction
            if (rxn->specLHS1 != nullptr) s = rxn->specLHS1; else s = rxn->specLHS2;
            addFirstOrderKineticTerm(s->index, s->index, -(rxn->rateConstantForwardNormalized));
            if (rxn->specRHS1 != nullptr) addFirstOrderKineticTerm(s->index, rxn->specRHS1->index, rxn->rateConstantForwardNormalized); //
            if (rxn->specRHS2 != nullptr) addFirstOrderKineticTerm(s->index, rxn->specRHS2->index, rxn->rateConstantForwardNormalized); //
        }
        // BACKWARD REACTIONS
        if (rxn->specRHS1 != nullptr && rxn->specRHS2 != nullptr)
        { // second-order backward reaction
            addSecondOrderKineticTerm(rxn->specRHS1->index, rxn->specRHS2->index, rxn->specRHS1->index, -(rxn->rateConstantBackwardNormalized));
            addSecondOrderKineticTerm(rxn->specRHS1->index, rxn->specRHS2->index, rxn->specRHS2->index, -(rxn->rateConstantBackwardNormalized));
            if (rxn->specLHS1 != nullptr) addSecondOrderKineticTerm(rxn->specRHS1->index, rxn->specRHS2->index, rxn->specLHS1->index, rxn->rateConstantBackwardNormalized);
            if (rxn->specLHS2 != nullptr) addSecondOrderKineticTerm(rxn->specRHS1->index, rxn->specRHS2->index, rxn->specLHS2->index, rxn->rateConstantBackwardNormalized);
        }
        else
        { // first-order backward reaction
            if (rxn->specRHS1 != nullptr) s = rxn->specRHS1; else s = rxn->specRHS2;
            addFirstOrderKineticTerm(s->index, s->index, -(rxn->rateConstantBackwardNormalized));
            if (rxn->specLHS1 != nullptr) addFirstOrderKineticTerm(s->index, rxn->specLHS1->index, rxn->rateConstantBackwardNormalized); //
            if (rxn->specLHS2 != nullptr) addFirstOrderKineticTerm(s->index, rxn->specLHS2->index, rxn->rateConstantBackwardNormalized); //
        }
    }
}

void Simulation::addFirstOrderKineticTerm(size_t spec1, size_t spec2, double rate)
{
    if (abs(rate) > MIN_RATE)
    {
        for(size_t x = 1; x < numGridPoints; x++)
        {
            addToMatrix(x+spec2*numGridPoints, x+spec1*numGridPoints, rate*deltaX*deltaX*paramGamma2i[x]);
        }
    }
}

void Simulation::addSecondOrderKineticTerm(size_t spec1, size_t spec2, size_t spec3, double rate)
{
    if (abs(rate) > MIN_RATE)
    {
        // spec1+spec2 react together, producing (at least) spec3 at rate rate
        secondOrderKinetics.emplace_back( make_tuple(spec1, spec2, spec3, rate) );
    }
}

void Simulation::addRedoxToMatrix()
{
    size_t oxidx, redidx, specidx;
    double p, Kred, Kox;

    for (auto redox: sys->vecRedox)
    {
        p = static_cast<double>(redox->numberElectrons) * f * (exper->initialPotential - redox->standardPotential); // normalized potential

        Kred = redox->rateConstantHeteroNormalized * exp(-redox->alpha * p); // K_red B-V kinetics ("Kf")
        Kox = redox->rateConstantHeteroNormalized * exp((1.0 - redox->alpha) * p); // K_ox B-V kinetics ("Kb")

        oxidx = redox->specOxidized->index*numGridPoints;
        redidx = redox->specReduced->index*numGridPoints;

        // add terms per redox reaction:
        addToMatrix(oxidx, redidx, - deltaX * Kox / redox->specOxidized->normalizedDiffusionConstant);  //B0
        addToMatrix(redidx, oxidx, - deltaX * Kred / redox->specReduced->normalizedDiffusionConstant); //B0
        // add to diagonal terms:
        speciesInRedox[redox->specOxidized->index] = true;
        matrixB0DiagonalTerms[redox->specOxidized->index] += deltaX * Kred / redox->specOxidized->normalizedDiffusionConstant;
        speciesInRedox[redox->specReduced->index] = true;
        matrixB0DiagonalTerms[redox->specReduced->index] += deltaX * Kox / redox->specReduced->normalizedDiffusionConstant;

        // set current contribution matrix:
        currentContributionMatrix.coeffRef(static_cast<Eigen::Index>(redox->specOxidized->index), static_cast<Eigen::Index>(redox->index)) = 1.0;
        currentContributionMatrix.coeffRef(static_cast<Eigen::Index>(redox->specReduced->index), static_cast<Eigen::Index>(redox->index)) = -1.0;
    }

    // add zero flux condition for all species NOT in redox step:
    for (auto spec: sys->vecSpecies)
    {
        specidx = spec->index*numGridPoints;

        if (!speciesInRedox[spec->index]) // species does not participate in any redox step
        {
            for (size_t x = 0; x < numCurrentPoints; x++)
            {
                // instead of zero flux condition function:
                addToMatrix(specidx, specidx+x, coeffBeta0[x]/deltaX);
            }
        }
        else // active redox species
        {
            addToMatrix(specidx, specidx, 1.0 + matrixB0DiagonalTerms[spec->index]);  // B0
            addToMatrix(specidx, specidx+1, -1.0);  // B1
        }
    }
}

void Simulation::addBICoeffsToMatrix()
{
    size_t row, relidx_max;
    for (auto spec: sys->vecSpecies)
    {
        for (size_t x = 1; x < numGridPoints; x++) // since x == 0 corresponds to the boundary condition
        {
            if (x < numGridPoints - numDerivPoints + 2)
                relidx_max = numDerivPoints-1;
            else
                relidx_max = numGridPoints-x;
            row = x + spec->index * numGridPoints; //row in matrixA
            for (size_t relidx = 0; relidx < relidx_max+1; relidx++)
            {
                addToMatrix(row, relidx+row-1, coeffMatrixN2(x, static_cast<int>(relidx)-1, spec->normalizedDiffusionConstant));
            }
        }
    }
}

double Simulation::coeffMatrixN2(size_t x, int relative_position, double normalizedDiffusionConstant)
{
    size_t coeffidx = static_cast<size_t>(relative_position+1); // relative_position >= -1
    double coeff = coeffAlpha[coeffidx] + el->electrodeGeometryFactor*deltaX*paramGammai[x]/(paramR0+gridCoordinate[x])*coeffBeta[coeffidx];
    coeff *= normalizedDiffusionConstant; // Compton page 88/89
    if (relative_position == 0) { coeff -= paramGamma2i[x]/paramLambda; }

    return coeff;
}

void Simulation::addToMatrix(size_t row, size_t col, double value)
{
    tripletContainerMatrixA.emplace_back( Triplet(static_cast<int>(row), static_cast<int>(col), value) );
}

void Simulation::createMatrix()
{
    // set initial values in matrix:
    matrixA.setFromTriplets(tripletContainerMatrixA.begin(), tripletContainerMatrixA.end());
}
