#include "simulationcore.h"

/*===============================================================================================
 * SIZING (PARAMETERS, COORDINATES, COEFFICIENTS)
 *=============================================================================================*/

double Sizing::initialize(System *sys, Electrode *el, Environment *env, Experiment *exper, ostream &output_stream)
{
    // some parameters
    f = CONST_F / (CONST_R * env->temperature);
    
    /* system dimensioning according to Compton (Understanding Voltammetry)
     * deltaTheta is a fixed value, taken from Compton
     */
    // calculate dimensionless parameters:
    const double sigma = exper->scanRate * f * std::pow(el->epsilon, 2) / sys->getDiffMax(); // scan rate
    /* kinetic parameter F:
     * to increase accuracy, we need to make deltaX smaller to accomodate for the fastest chemical reaction
     * the factor F in deltaX = pow(10.0, -F-0.5*log10(sigma)); should e.g. go from 2.2 to 4.2
     */
    double F;
    const double logRate = log10(sys->getRateMax());
    if (logRate < minLogRate) { F = minF; }
    else if (logRate <= maxLogRate) { F = minF + (logRate - minLogRate) * (maxF - minF) / (maxLogRate - minLogRate); }
    else { F = maxF; }
    // distance and time step:
    deltaX = pow(10.0, -F) / sqrt(sigma); // Compton (page 64) leads to: dX = 10^(-2.2-0.5*sigma), deltaX = size of smallest grid element [1]
    deltaT = deltaTheta / sigma; // deltaT = size of smallest time step [1]
    // electrode and current factors & parameters:
    currentFromFlux = el->epsilon * CONST_F * sys->getDiffMax() * sys->getConcMax();
    
    // create expanding grid:
    const double maxX = 6.0*sqrt(exper->totalTime() * f * exper->scanRate / sigma); // maxX = maximum distance from electrode [1]
    numGridPoints = static_cast<int>(std::ceil( std::log((maxX/deltaX * (paramGamma-1.0))+1.0) / std::log(paramGamma) )) + 1;
    // set system dimensions:
    numCurrentPoints = min(numCurrentPoints, numGridPoints); // truncate number of current points if it is larger than number of grid points
    numSpecies = sys->vecSpecies.size();
    numRedox = sys->vecRedox.size();
    numReactions = sys->vecReactions.size();
    mainDimension = numSpecies*numGridPoints;
    
    // initialize distance in grid, expansion factor gamma and concentration vectors:
    paramGamma2i2dX.resize(numGridPoints);
    for (size_t x = 0; x < numGridPoints; x++)
        paramGamma2i2dX[x] = std::pow(paramGamma, 2*x) * std::pow(deltaX, 2);
    
    // obtain flux condition coefficients:
    coeffBeta0.resize(numCurrentPoints);
    for (size_t i = 0; i < numCurrentPoints; i++)
        coeffBeta0[i] = Beta_N_1(numCurrentPoints, i, paramGamma);
    
    // obtain (N,2) coefficients "a" and populate matrices (one for each species)
    // see also Molina (10.1016/j.cplett.2015.11.011) equation (14)
    coeffN2.resize(numSpecies);
    double coeff, gridCoordinate = 0, gammaxdx, alpha, beta;
    const double paramR0 = el->epsilon / sqrt(sys->getDiffMax() * exper->totalTime());
    for (size_t s = 0; s < numSpecies; s++)
    {
        coeffN2[s] = Eigen::MatrixXd(numGridPoints, numDerivPoints);
        for (size_t x = 0; x < numGridPoints; x++)
        {
            gammaxdx = std::pow(paramGamma, x) * deltaX;
            for (size_t d = 0; d < numDerivPoints; d++)
            {
                alpha = Alpha_N_2(numDerivPoints, d-1, paramGamma);
                beta  =  Beta_N_2(numDerivPoints, d-1, paramGamma);
                
                coeff = alpha + beta*el->electrodeGeometryFactor*gammaxdx/(paramR0+gridCoordinate);
                coeff *= sys->vecSpecies[s]->getDiffNorm(); // Compton page 88/89
                coeff -= (d == 1) ? paramGamma2i2dX[x]/deltaT : 0.0;
                coeffN2[s](x, d) = coeff;
            }
            gridCoordinate += gammaxdx;
        }
    }
    
    // print scaling parameters:
    /*
    output_stream << "Initial concentrations after equilibration: (";
    for (auto spec: sys->vecSpecies)
        output_stream << "[" << spec->getName() << "] = " << spec->getConcNormEquil() << ((spec!=sys->vecSpecies.back())?", ":"");
    output_stream << ")<br />" << endl;
    output_stream << "V<sub>step</sub> [V] = " << deltaTheta/f << ", t<sub>max</sub> [s] = " << totalTime << ", " << endl;
    output_stream << "r<sub>max</sub> [1/s] = " << sys->getRateMax() << ", F = " << F << "<br />" << endl;
    output_stream << "&sigma; = " << sigma << ", &theta;<sub>max</sub> = " << totalTheta << ", &Delta;&theta; = " << deltaTheta << ", " << endl;
    output_stream << "X<sub>max</sub> = " << maxX << ", &Delta;X = " << deltaX << ", T<sub>max</sub> = " << maxT << ", &Delta;T = " << deltaT << "<br />" << endl;

    // output parametrisation of grid:
    output_stream << "Number of grid points = " << numGridPoints << "<br />" << endl;
    output_stream << "Number of current points = " << numCurrentPoints << endl;
    output_stream << "Number of deriv points = " << numDerivPoints << endl;
    */
    
    return deltaTheta/f;
}

/*===============================================================================================
 * CORE
 *=============================================================================================*/

void Core::initialize(Sizing *new_sz, System *new_sys, double ip)
{
    sz = new_sz;
    sys = new_sys;
    
    // initialize matrix system:
    msys.initialize(sz);
    // initialize supporting matrices and vectors:
    initVectorsEtc();
    // fill matrix with terms and coefficients:
    addRedoxToMatrix(ip);
    addBICoeffsToMatrix();
    addKineticsToMatrix();
    // create matrix:
    msys.createMatrix();
}

void Core::solveSystem(double potential)
{
    gridConcentrationPrevious = gridConcentration; // store concentrations
    
    updateIndependentTerms(); // update independent terms in matrix
    updateRedoxInMatrix(potential); // update redox terms in matrix
    updateKineticsInMatrix(1.0); // add 2nd order kinetic terms (Laasonen) in matrix
    
    msys.solve(&independentTerms, &gridConcentration); // solve the matrix system
    
    updateKineticsInMatrix(-1.0); // subtract 2nd order kinetic terms (Laasonen) from matrix
}

/////////////////////////////////////////////

void Core::updateKineticsInMatrix(double factor)
{
    double rate;
    size_t spec1, spec2, spec3;

    for (auto kinTerm: secondOrderKinetics)
    {
        spec1 = get<0>(kinTerm);
        spec2 = get<1>(kinTerm);
        spec3 = get<2>(kinTerm);

        for (size_t x = 1; x < sz->numGridPoints; x++)
        {
            rate = factor * get<3>(kinTerm) * sz->paramGamma2i2dX[x];
            
            msys.addToValue( spec3, x, spec1, x, rate * gridConcentration(x, spec2) );
            msys.addToValue( spec3, x, spec2, x, rate * gridConcentration(x, spec1) );
            
            independentTerms(x, spec3) += rate * gridConcentration(x, spec1) * gridConcentration(x, spec2);
        }
    }
}

void Core::updateIndependentTerms()
{
    for (size_t s = 0; s < sz->numSpecies; s++)
    {
        //Zero at matrix rows storing surface conditions
        independentTerms(0, s) = 0.0;

        for (int x = 1; x < sz->numGridPoints; x++)
        {
            /* independentTerms is the RHS ('b') of the matrix equation A*x=b, where 'x' is the concentration at T+dT (the 'next' concentration)
             * the system of equations is ('a' are coefficients, 'x(i,s)' the next concentrations at grid point i for species s):
             * a(-1)*x(i-1,s) + a(0)*x(i,s) + a(+1)*x(i+1,s) + a(+2)*x(i+2,s) + (etc.) = b(i,s)
             */
            independentTerms(x, s) = -gridConcentration(x, s) * sz->paramGamma2i2dX[x]/sz->deltaT; // see also Molina (10.1016/j.cplett.2015.11.011) equation (14)
            
            // however, when i>numGridPoints, then x(i+i) reduces to the (never changing) bulk concentration and we add them into 'b', because they are now known:
            if (x > sz->numGridPoints - sz->numDerivPoints + 1)
                for (int d = sz->numGridPoints-x+1; d < sz->numDerivPoints; d++)
                    independentTerms(x, s) -= sz->coeffN2[s](x, d) * sys->vecSpecies[s]->getConcNormEquil();
        }
    }
}

void Core::updateRedoxInMatrix(double potential)
{
    double p, Kred, Kox;

    // reset diagonal terms:
    fill(matrixB0DiagonalTerms.begin(), matrixB0DiagonalTerms.end(), 0.0);

    for (auto redox: sys->vecRedox)
    {
        p = static_cast<double>(redox->getNe()) * (sz->f) * (potential - redox->getE0()); // normalized potential

        Kred = (sz->deltaX) * redox->getKeNorm() * exp(-redox->getAlpha() * p); // B-V kinetics
        Kox =  (sz->deltaX) * redox->getKeNorm() * exp((1.0 - redox->getAlpha()) * p); // B-V kinetics

        // add terms to matrix per redox reaction:
        msys.setValue(  redox->getSpecOx()->getIndex(), 0, redox->getSpecRed()->getIndex(), 0, -Kox / redox->getSpecOx()->getDiffNorm() );   // B0
        msys.setValue( redox->getSpecRed()->getIndex(), 0, redox->getSpecOx()->getIndex(),  0, -Kred / redox->getSpecRed()->getDiffNorm() ); // B0

        // add to diagonal terms:
        matrixB0DiagonalTerms[redox->getSpecOx()->getIndex()] += Kred / redox->getSpecOx()->getDiffNorm();
        matrixB0DiagonalTerms[redox->getSpecRed()->getIndex()] += Kox / redox->getSpecRed()->getDiffNorm();
    }

    // add diagonal terms to matrix:
    for (size_t s = 0; s < sz->numSpecies; s++)
    {
        if (speciesInRedox[s]) // active redox species
        {
            msys.setValue( s, 0, s, 0,  1.0 + matrixB0DiagonalTerms[s] ); // B0
            msys.setValue( s, 0, s, 1, -1.0 ); // B1
        }
    }
}

double Core::calcCurrentFromFlux()
{
    double totalflux, speciesflux;
    
    // determine flux per species:
    for (auto spec: sys->vecSpecies)
    {
        speciesflux = 0;
        for(size_t x = 0; x < sz->numCurrentPoints; x++)
            speciesflux += sz->coeffBeta0[x] * gridConcentration(x, spec->getIndex());
        currentContributionSpeciesFlux[spec->getIndex()] = speciesflux * spec->getDiffNorm();
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
        totalflux -= currentContributionRedoxFlux[redox->getIndex()] * static_cast<double>(redox->getNe());
    totalflux *= (sz->currentFromFlux) / (sz->deltaX); // totalflux := dC/dX

    // calculate current from flux:
    return totalflux;
}

/////////////////////////////////////////////

void Core::initVectorsEtc()
{
    // initialize independent terms and grid concentrations:
    independentTerms = Eigen::MatrixXd::Zero(sz->numGridPoints, sz->numSpecies);
    gridConcentration = Eigen::MatrixXd(sz->numGridPoints, sz->numSpecies);
    gridConcentrationPrevious = Eigen::MatrixXd::Zero(sz->numGridPoints, sz->numSpecies);
    // load in starting concentrations:
    for (size_t x = 0; x < sz->numGridPoints; x++)
        for (auto spec: sys->vecSpecies)
            gridConcentration(x, spec->getIndex()) = spec->getConcNormEquil();

    // initilize redox/current vectors & matrix:
    speciesInRedox.clear();
    speciesInRedox.resize(sz->numSpecies, false);
    matrixB0DiagonalTerms.clear();
    matrixB0DiagonalTerms.resize(sz->numSpecies, 0.0);
    currentContributionMatrix = Eigen::MatrixXd::Zero(sz->numSpecies, sz->numSpecies);
    currentContributionSpeciesFlux = Eigen::VectorXd::Zero(sz->numSpecies);
    currentContributionRedoxFlux = Eigen::VectorXd::Zero(sz->numRedox);
}

void Core::addKineticsToMatrix()
{
    for (auto rxn: sys->vecReactions) {
        addHalfReactionKineticsToMatrix(rxn->getSpecLHS1(), rxn->getSpecLHS2(),
                                        rxn->getSpecRHS1(), rxn->getSpecRHS2(), rxn->getKfNorm()); // FORWARD REACTIONS
        addHalfReactionKineticsToMatrix(rxn->getSpecRHS1(), rxn->getSpecRHS2(),
                                        rxn->getSpecLHS1(), rxn->getSpecLHS2(), rxn->getKbNorm()); // BACKWARD REACTIONS
    }
}

void Core::addHalfReactionKineticsToMatrix(Species *f1, Species *f2, Species *b1, Species *b2, double normrate)
{
    if (f1 != nullptr && f2 != nullptr)
    {
        // second-order reaction
        addSecondOrderKineticTerm(f1, f2, f1, -normrate);
        addSecondOrderKineticTerm(f1, f2, f2, -normrate);
        addSecondOrderKineticTerm(f1, f2, b1, normrate);
        addSecondOrderKineticTerm(f1, f2, b2, normrate);
    }
    else if (f1 != nullptr || f2 != nullptr)
    {
        // first-order reaction
        Species *f = (f1 != nullptr) ? f1 : f2;
        addFirstOrderKineticTerm(f, f, -normrate);
        addFirstOrderKineticTerm(f, b1, normrate);
        addFirstOrderKineticTerm(f, b2, normrate);
    }
    else
    {
        // reaction has either no products or no reactants; is that a problem? throw error?
    }
}

void Core::addFirstOrderKineticTerm(Species *spec1, Species *spec2, double normrate)
{
    if (abs(normrate) > MIN_RATE && spec2 != nullptr)
        for(size_t x = 1; x < sz->numGridPoints; x++)
            msys.setValue(spec2->getIndex(), x, spec1->getIndex(), x, normrate*sz->paramGamma2i2dX[x]);
}

void Core::addSecondOrderKineticTerm(Species *spec1, Species *spec2, Species *spec3, double normrate)
{
    if (abs(normrate) > MIN_RATE && spec3 != nullptr)
        secondOrderKinetics.push_back( make_tuple(spec1->getIndex(), spec2->getIndex(), spec3->getIndex(), normrate) );
}

void Core::addRedoxToMatrix(double ip)
{
    for (auto redox: sys->vecRedox)
    {
        // keep track of which species partake in redox steps:
        speciesInRedox[redox->getSpecOx()->getIndex()] = true;
        speciesInRedox[redox->getSpecRed()->getIndex()] = true;
    
        // set current contribution matrix:
        currentContributionMatrix.coeffRef(redox->getSpecOx()->getIndex(),  redox->getIndex()) = 1.0;
        currentContributionMatrix.coeffRef(redox->getSpecRed()->getIndex(), redox->getIndex()) = -1.0;
    }
    
    // instead of zero flux condition function:
    for (size_t s = 0; s < sz->numSpecies; s++)
        if (!speciesInRedox[s])
            for (size_t x = 0; x < sz->numCurrentPoints; x++)
                msys.setValue(s, 0, s, x, sz->coeffBeta0[x]/sz->deltaX);

    updateRedoxInMatrix(ip);
}

void Core::addBICoeffsToMatrix() // Backwards Implicit coefficients
{
    size_t relidx_max;
    for (size_t s = 0; s < sz->numSpecies; s++)
    {
        for (size_t x = 1; x < sz->numGridPoints; x++) // since x == 0 corresponds to the boundary condition
        {
            // the "relative index" goes from -1 to numDerivPoints-1, but is bounded by the size of the grid:
            relidx_max = std::min(sz->numDerivPoints-1, sz->numGridPoints-x);
            
            for (int relidx = -1; relidx < static_cast<int>(relidx_max); relidx++)
                msys.setValue(s, x, s, x+relidx, sz->coeffN2[s](x, relidx+1));
        }
    }
}

/*===============================================================================================
 * MATRIXSYSTEM
 *=============================================================================================*/

using namespace std::placeholders;

void MatrixSystem::initialize(Sizing *new_sz)
{
    sz = new_sz;
    const size_t numNonZeroesTotal = 15*sz->mainDimension;  // estimation (15 > 2*numDerivPoints)
    
    vecb = Eigen::VectorXd::Zero(sz->mainDimension);
    vecx = Eigen::VectorXd::Zero(sz->mainDimension);
    
    // initialize Triplet container
    tripletContainerMatrixA.clear();
    tripletContainerMatrixA.reserve(numNonZeroesTotal);

    // clear and resize matrix:
    matrixA.setZero(); // removes all non-zeros
    matrixA.resize(0,0);
    matrixA.data().squeeze(); // removes all items from matrix
    matrixA.resize(sz->mainDimension, sz->mainDimension);
    matrixA.reserve(numNonZeroesTotal);
    matrixPatternAnalyzed = false;
    
    // change setValue function pointer:
    setValue = std::bind(&MatrixSystem::setValueTriplet, this, _1, _2, _3, _4, _5);
}

void MatrixSystem::solve(Eigen::MatrixXd *independentTerms, Eigen::MatrixXd *gridConcentration)
{
    // pattern needs to be analyzed once after matrix changed:
    if (!matrixPatternAnalyzed) { sparseMatrixSolver.analyzePattern(matrixA); matrixPatternAnalyzed = true; }

    // flatten matrix to vector (because solver uses 2D*1D=1D sizing):
    vecb = Eigen::Map<Eigen::VectorXd>(independentTerms->data(), sz->mainDimension);
    // this is where the magic happens:
    sparseMatrixSolver.factorize(matrixA);
    vecx = sparseMatrixSolver.solve(vecb);
    // reshape vector into matrix:
    *gridConcentration = Eigen::Map<Eigen::MatrixXd>(vecx.data(), sz->numGridPoints, sz->numSpecies);
}

void MatrixSystem::createMatrix()
{
    // create matrix from triplets:
    matrixA.setFromTriplets(tripletContainerMatrixA.begin(), tripletContainerMatrixA.end());
    
    // change setValue function pointer:
    setValue = std::bind(&MatrixSystem::setValueMatrix, this, _1, _2, _3, _4, _5);
}
