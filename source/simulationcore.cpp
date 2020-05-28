#include "simulationcore.h"

/*===============================================================================================
 * SIZING
 *=============================================================================================*/

double Sizing::initialize(System *sys, Electrode *el, Environment *env, Experiment *exper, ostream &output_stream)
{
    // some parameters
    f = CONST_F / (CONST_R * env->temperature);
    double totalTime = exper->totalTime();
    double totalTheta = totalTime * f * exper->scanRate;
    
    /* system dimensioning according to Compton (Understanding Voltammetry)
     * deltaTheta is a fixed value, taken from Compton
     * to increase accuracy, we need to decrease deltaX, which doesn't need a lot of extra time to solve!!!
     * the factor F in deltaX = pow(10.0, -F-0.5*log10(sigma)); should e.g. go from 2.2 to 4.2
     */
    double deltaE = deltaTheta/f;
    double F, logRate = log10(sys->getRateMax());
    if (logRate < minLogRate) { F = minF; }
    else if (logRate <= maxLogRate) { F = minF + (logRate - minLogRate) * (maxF - minF) / (maxLogRate - minLogRate); }
    else { F = maxF; }
    double sigma = el->epsilon * el->epsilon / sys->getDiffMax() * f * exper->scanRate; // dimensionless scan rate
    double maxT = totalTheta / sigma; // maxT = total experiment time [1]
    double  deltaT = deltaTheta / sigma; // deltaT = size of smallest time step [1]
    double maxX = 6.0*sqrt(maxT); // maxX = maximum distance from electrode [1]
    deltaX = pow(10.0, -F) / sqrt(sigma); // Compton (page 64) leads to: dX = 10^(-2.2-0.5*sigma), deltaX = size of smallest grid element [1]
    paramLambda = deltaT / (deltaX * deltaX);
    currentFromFlux = el->epsilon * CONST_F * sys->getDiffMax() * sys->getConcMax();
    electrodeGeometryFactor = el->electrodeGeometryFactor;
    paramR0 = el->epsilon / sqrt(sys->getDiffMax() * totalTime);

    // create expanding grid and set dimensions:
    numGridPoints = 1;
    do { numGridPoints++; } while (deltaX < maxX * (paramGamma - 1.0) / ( pow(paramGamma, numGridPoints-1) - 1.0 ));
    numCurrentPoints = min(numCurrentPoints, numGridPoints); // truncate number of current points if it is larger than number of grid points
    numSpecies = sys->vecSpecies.size();
    numRedox = sys->vecRedox.size();
    numReactions = sys->vecReactions.size();
    
        // print scaling parameters:
    //output_stream << "Diff coeff [m2/s]: max = " << sys.getDiffMax() << endl;
    //output_stream << "Conc [mol/m3]: max = " << sys.getConcMax() << endl;
    //output_stream << "Electrode [m]: epsilon = " << el.epsilon << ", area = " << el.electrodeArea << endl;
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
    
    return deltaE;
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
            rate = factor * get<3>(kinTerm) * (sz->deltaX)*(sz->deltaX)*paramGamma2i[x];
            
            msys.addToValue( spec3, x, spec1, x, rate * gridConcentration(spec2, x) );
            msys.addToValue( spec3, x, spec2, x, rate * gridConcentration(spec1, x) );
            
            independentTerms(spec3, x) += rate * gridConcentration(spec1, x) * gridConcentration(spec2, x);
        }
    }
}

void Core::updateIndependentTerms()
{
    size_t EndNormal = sz->numGridPoints - sz->numDerivPoints + 2;

    for (auto spec: sys->vecSpecies)
    {
        //Zero at matrix rows storing surface conditions
        independentTerms(spec->getIndex(), 0) = 0.0;

        for (size_t x = 1; x < sz->numGridPoints; x++)
        {
            /* independentTerms is the RHS ('b') of the matrix equation A*x=b, where 'x' is the concentration at T+dT (the 'next' concentration)
             * the system of equations is ('a' are coefficients, 'x(i,s)' the next concentrations at grid point i for species s):
             * a(-1)*x(i-1,s) + a(0)*x(i,s) + a(+1)*x(i+1,s) + a(+2)*x(i+2,s) + (etc.) = b(i,s)
             */
            independentTerms(spec->getIndex(), x) = -gridConcentration(spec->getIndex(), x) * paramGamma2i[x]/(sz->paramLambda);
            // however, when i>numGridPoints, then x(i+i) reduces to the (never changing) bulk concentration and we add them into 'b', because they are now known:
            if (x >= EndNormal)
            {
                for (int jj = 0; jj <= static_cast<int>(x - EndNormal); jj++)
                {
                    independentTerms(spec->getIndex(), x) -= coeffMatrixN2(x, static_cast<int>(sz->numDerivPoints)-jj-2, spec->getDiffNorm()) * spec->getConcNormEquil();
                }
            }
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
        msys.changeValue( redox->getSpecOx()->getIndex(), 0,
                          redox->getSpecRed()->getIndex(), 0,
                          -Kox / redox->getSpecOx()->getDiffNorm() );   // B0
        msys.changeValue( redox->getSpecRed()->getIndex(), 0,
                          redox->getSpecOx()->getIndex(), 0, 
                          -Kred / redox->getSpecRed()->getDiffNorm() ); // B0

        // add to diagonal terms:
        matrixB0DiagonalTerms[redox->getSpecOx()->getIndex()] += Kred / redox->getSpecOx()->getDiffNorm();
        matrixB0DiagonalTerms[redox->getSpecRed()->getIndex()] += Kox / redox->getSpecRed()->getDiffNorm();
    }

    // add diagonal terms to matrix:
    for (size_t s = 0; s < sz->numSpecies; s++)
    {
        if (speciesInRedox[s]) // active redox species
        {
            msys.changeValue( s, 0, s, 0,  1.0 + matrixB0DiagonalTerms[s] ); // B0
            msys.changeValue( s, 0, s, 1, -1.0 ); // B1
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
        for(size_t kk = 0; kk < sz->numCurrentPoints; kk++)
        {
            speciesflux += coeffBeta0[kk] * gridConcentration(spec->getIndex(), kk);
        }
        currentContributionSpeciesFlux[static_cast<Eigen::Index>(spec->getIndex())] = speciesflux * spec->getDiffNorm();
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
        totalflux -= currentContributionRedoxFlux[static_cast<Eigen::Index>(redox->getIndex())]
                     * static_cast<double>(redox->getNe());
    }
    totalflux /= (sz->deltaX); // totalflux := dC/dX

    // calculate current from flux:
    return totalflux * (sz->currentFromFlux);
}

/////////////////////////////////////////////

void Core::initVectorsEtc()
{
    // initialize independent terms and b & x vectors:
    independentTerms = Eigen::MatrixXd::Zero(sz->numSpecies, sz->numGridPoints);
    gridConcentration = Eigen::MatrixXd(sz->numSpecies, sz->numGridPoints);
    gridConcentrationPrevious = Eigen::MatrixXd::Zero(sz->numSpecies, sz->numGridPoints);
    
    // initialize distance in grid, expansion factor gamma and concentration vectors:
    gridCoordinate.resize(sz->numGridPoints);
    paramGammai.resize(sz->numGridPoints);
    paramGamma2i.resize(sz->numGridPoints);
    for (size_t x = 0; x < sz->numGridPoints; x++)
    {
        paramGammai[x] = pow(sz->paramGamma, x);
        paramGamma2i[x] = paramGammai[x]*paramGammai[x];
        gridCoordinate[x] = (x == 0) ? 0.0 : gridCoordinate[x-1] + (sz->deltaX)*paramGammai[x-1];

        for (auto spec: sys->vecSpecies)
            gridConcentration(spec->getIndex(), x) = spec->getConcNormEquil();
    }
    
    // obtain flux condition coefficients:
    coeffBeta0.resize(sz->numCurrentPoints);
    for(size_t i = 0; i < sz->numCurrentPoints; i++)
        coeffBeta0[i] = Beta_N_1(static_cast<int>(sz->numCurrentPoints), static_cast<int>(i), sz->paramGamma);
    
    // obtain diffusion coefficients:
    coeffAlpha.resize(sz->numDerivPoints);
    coeffBeta.resize(sz->numDerivPoints);
    for(size_t d = 0; d < sz->numDerivPoints; d++)
    {
        coeffAlpha[d] = Alpha_N_2(static_cast<int>(sz->numDerivPoints), static_cast<int>(d)-1, sz->paramGamma);
        coeffBeta[d]  =  Beta_N_2(static_cast<int>(sz->numDerivPoints), static_cast<int>(d)-1, sz->paramGamma);
    }

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
        // reaction has no products or reactants; is that a problem?
    }
}

void Core::addFirstOrderKineticTerm(Species *spec1, Species *spec2, double normrate)
{
    if (abs(normrate) > MIN_RATE && spec2 != nullptr)
        for(size_t x = 1; x < sz->numGridPoints; x++)
            msys.addToMatrix(spec2->getIndex(), x, spec1->getIndex(), x, normrate*(sz->deltaX)*(sz->deltaX)*paramGamma2i[x]);
}

void Core::addSecondOrderKineticTerm(Species *spec1, Species *spec2, Species *spec3, double normrate)
{
    if (abs(normrate) > MIN_RATE && spec3 != nullptr)
        secondOrderKinetics.push_back( make_tuple(spec1->getIndex(), spec2->getIndex(), spec3->getIndex(), normrate) );
}

void Core::addRedoxToMatrix(double ip)
{
    double p, Kred, Kox;

    for (auto redox: sys->vecRedox)
    {
        p = static_cast<double>(redox->getNe()) * (sz->f) * (ip - redox->getE0()); // normalized potential

        Kred = (sz->deltaX) * redox->getKeNorm() * exp(      -redox->getAlpha()  * p); // K_red B-V kinetics ("Kf")
        Kox  = (sz->deltaX) * redox->getKeNorm() * exp((1.0 - redox->getAlpha()) * p); // K_ox  B-V kinetics ("Kb")

        // add terms per redox reaction:
        msys.addToMatrix(redox->getSpecOx()->getIndex(),  0,
                         redox->getSpecRed()->getIndex(), 0,
                         -Kox / redox->getSpecOx()->getDiffNorm());  //B0
        msys.addToMatrix(redox->getSpecRed()->getIndex(), 0,
                         redox->getSpecOx()->getIndex(),  0,
                         -Kred / redox->getSpecRed()->getDiffNorm()); //B0
        
        // keep track of which species partake in redox steps:
        speciesInRedox[redox->getSpecOx()->getIndex()] = true;
        speciesInRedox[redox->getSpecRed()->getIndex()] = true;
        
        // add to diagonal terms:
        matrixB0DiagonalTerms[redox->getSpecOx()->getIndex()] += Kred / redox->getSpecOx()->getDiffNorm();
        matrixB0DiagonalTerms[redox->getSpecRed()->getIndex()] += Kox / redox->getSpecRed()->getDiffNorm();

        // set current contribution matrix:
        currentContributionMatrix.coeffRef(static_cast<Eigen::Index>(redox->getSpecOx()->getIndex()), static_cast<Eigen::Index>(redox->getIndex())) = 1.0;
        currentContributionMatrix.coeffRef(static_cast<Eigen::Index>(redox->getSpecRed()->getIndex()), static_cast<Eigen::Index>(redox->getIndex())) = -1.0;
    }

    // add zero flux condition for all species NOT in redox step:
    for (size_t s = 0; s < sz->numSpecies; s++)
    {
        if (!speciesInRedox[s]) // species does not participate in any redox step
        {
            for (size_t x = 0; x < sz->numCurrentPoints; x++)
            {
                // instead of zero flux condition function:
                msys.addToMatrix(s, 0, s, x, coeffBeta0[x]/(sz->deltaX));
            }
        }
        else // active redox species
        {
            msys.addToMatrix(s, 0, s, 0,  1.0 + matrixB0DiagonalTerms[s]);  // B0
            msys.addToMatrix(s, 0, s, 1, -1.0);  // B1
        }
    }
}

void Core::addBICoeffsToMatrix() // Backwards Implicit coefficients
{
    size_t relidx_max;
    for (auto spec: sys->vecSpecies)
    {
        for (size_t x = 1; x < sz->numGridPoints; x++) // since x == 0 corresponds to the boundary condition
        {
            if (x < sz->numGridPoints - sz->numDerivPoints + 2)
                relidx_max = sz->numDerivPoints-1;
            else
                relidx_max = sz->numGridPoints-x;
            
            for (size_t relidx = 0; relidx < relidx_max+1; relidx++)
            {
                msys.addToMatrix(spec->getIndex(), x, spec->getIndex(), x+relidx-1,
                                 coeffMatrixN2(x, static_cast<int>(relidx)-1, spec->getDiffNorm()));
            }
        }
    }
}

double Core::coeffMatrixN2(size_t x, int relative_position, double diffnorm)
{
    size_t coeffidx = static_cast<size_t>(relative_position+1); // relative_position >= -1
    
    double coeff = coeffAlpha[coeffidx] + (sz->electrodeGeometryFactor)*(sz->deltaX)*paramGammai[x]/((sz->paramR0)+gridCoordinate[x])*coeffBeta[coeffidx];
    coeff *= diffnorm; // Compton page 88/89
    coeff -= (relative_position == 0) ? paramGamma2i[x]/(sz->paramLambda) : 0.0;
    
    return coeff;
}

/*===============================================================================================
 * MATRIXSYSTEM
 *=============================================================================================*/

void MatrixSystem::initialize(Sizing *new_sz)
{
    sz = new_sz;
    size_t numNonZeroesTotal = 15*sz->numSpecies*sz->numGridPoints;  // estimation (15 > 2*numDerivPoints)
    
    vecb = Eigen::VectorXd::Zero(sz->numSpecies*sz->numGridPoints);
    vecx = Eigen::VectorXd::Zero(sz->numSpecies*sz->numGridPoints);
    
    // initialize Triplet container
    tripletContainerMatrixA.clear();
    tripletContainerMatrixA.reserve(numNonZeroesTotal);

    // clear and resize matrix:
    matrixA.setZero(); // removes all non-zeros
    matrixA.resize(0,0);
    matrixA.data().squeeze(); // removes all items from matrix
    matrixA.resize(static_cast<Eigen::Index>(sz->numSpecies*sz->numGridPoints),
                   static_cast<Eigen::Index>(sz->numSpecies*sz->numGridPoints));
    matrixA.reserve(static_cast<Eigen::Index>(numNonZeroesTotal));
    
    matrixPatternAnalyzed = false;
}

void MatrixSystem::solve(Eigen::MatrixXd *independentTerms, Eigen::MatrixXd *gridConcentration)
{
    // pattern needs to be analyzed once after matrix changed:
    if (!matrixPatternAnalyzed) { sparseMatrixSolver.analyzePattern(matrixA); matrixPatternAnalyzed = true; }
    
    sparseMatrixSolver.factorize(matrixA);

    // flatten matrix to vector:
    for (size_t s = 0; s < sz->numSpecies; s++)
        for (size_t x = 0; x < sz->numGridPoints; x++)
            vecb[static_cast<Eigen::Index>( s*(sz->numGridPoints)+x )] = (*independentTerms)(s, x);
    
    // this is where the magic happens:
    vecx = sparseMatrixSolver.solve(vecb);
    
    // reshape vector into matrix:
    for (size_t s = 0; s < sz->numSpecies; s++)
        for (size_t x = 0; x < sz->numGridPoints; x++)
            (*gridConcentration)(s, x) = vecx[static_cast<Eigen::Index>( s*(sz->numGridPoints)+x )];
}

void MatrixSystem::addToMatrix(size_t s1, size_t g1, size_t s2, size_t g2, double value)
{
    tripletContainerMatrixA.push_back( Triplet(static_cast<int>(s1*(sz->numGridPoints)+g1),
                                                  static_cast<int>(s2*(sz->numGridPoints)+g2), value) );
}

void MatrixSystem::createMatrix()
{
    // set initial values in matrix:
    matrixA.setFromTriplets(tripletContainerMatrixA.begin(), tripletContainerMatrixA.end());
}

void MatrixSystem::changeValue(size_t s1, size_t g1, size_t s2, size_t g2, double value)
{
    matrixA.coeffRef(static_cast<int>(s1*(sz->numGridPoints)+g1),
                     static_cast<int>(s2*(sz->numGridPoints)+g2)) = value;
}

void MatrixSystem::addToValue(size_t s1, size_t g1, size_t s2, size_t g2, double value)
{
    matrixA.coeffRef(static_cast<int>(s1*(sz->numGridPoints)+g1),
                     static_cast<int>(s2*(sz->numGridPoints)+g2)) += value;
}
 
