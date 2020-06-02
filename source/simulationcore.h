#ifndef SIMULATIONCORE_H
#define SIMULATIONCORE_H

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm> // integer max() & min()

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/Core>

const double MIN_RATE = 1.0e-6; // minimal rate > 1e-6 /s seems reasonable (as in: no rate will ever be smaller than this)
const double MIN_CONC = 1.0e-9; // minimal (initial!) concentration > 1 pM (because /1000) seems reasonable
const double MIN_AREA = 1.0e-18; // minimal electrode area > 1 nm2 seems reasonable

const double CONST_PI = 3.14159265358979323846;
const double CONST_F = 96485.33212;
const double CONST_R = 8.31446261815324;

#include "coefs_alpha_beta.h"
#include "electrodes.h"
#include "system.h"
#include "environment.h"
#include "experiment.h"

class Sizing
{
public:
    /*
     * PARAMETERS
     */
    // physical:
    double f; // F/RT
    // dimensionless parameters:
    double deltaTheta = 0.2; // potential step
    double paramGamma = 1.67; // grid expansion factor
    double minF = 2.2, maxF = 6.2, minLogRate = 3.0, maxLogRate = 7.0; // parameters that relate to homogeneous reaction rates
    double deltaX, deltaT; // smallest dimensionless distance & time
    double paramR0, electrodeGeometryFactor; // electrode parameters
    // conversion:
    double currentFromFlux; // to calculate current from species' flux
    // sizes of system parts:
    std::size_t numSpecies, numRedox, numReactions;
    std::size_t numGridPoints, numCurrentPoints = 5, numDerivPoints = 6; // number of points in grid, number of points for calculation of current, number of points for diffusion
    std::size_t mainDimension; // used heavily: numSpecies*numGridPoints
    /*
     * COEFFICIENTS & COORDINATES
     */
    std::vector<double> paramGamma2i2dX; // (grid expansion factor^x * deltaX)^2, pre-calculated for all points in grid
    std::vector<double> coeffBeta0; // coefficients for current (beta0)
    std::vector<Eigen::MatrixXd> coeffN2; // matrix storing all the pre-calculated diffusion coefficients
    
    /*
     * METHODS
     */
    double initialize(System*, Electrode*, Environment*, Experiment*, std::ostream &);
};

using Triplet = Eigen::Triplet<double>; //triplets (row, column, value) to fill sparse matrix
using Solver  = Eigen::SparseLU<Eigen::SparseMatrix<double>>;  //Turns out that SparseLU is actually fast enough, but only when running as Release :)

class MatrixSystem
{
private:
    Sizing* sz;
    
    std::vector<Triplet> tripletContainerMatrixA; // this container holds (row, col, value) triplets before matrixA is assembled
    Eigen::SparseMatrix<double> matrixA; // matrix to be inverted: Ax=b corresponds to matrixA * vecx = vecb
    Eigen::VectorXd vecb, vecx; // vectors, defininition see above
    Solver sparseMatrixSolver; // solver that inverts the sparse matrixA using LU factorization
    
    bool matrixPatternAnalyzed;
    
    inline std::size_t idx(std::size_t s, std::size_t x) { return (s*(sz->numGridPoints)+x); };
    inline void setValueMatrix(std::size_t s1, std::size_t g1, std::size_t s2, std::size_t g2, double value)
        { matrixA.coeffRef(idx(s1,g1), idx(s2,g2)) = value; };
    inline void setValueTriplet(std::size_t s1, std::size_t g1, std::size_t s2, std::size_t g2, double value)
        { tripletContainerMatrixA.push_back(Triplet( static_cast<int>(idx(s1,g1)), static_cast<int>(idx(s2,g2)), value)); };
public:
    void initialize(Sizing*);
    
    void solve(Eigen::MatrixXd*, Eigen::MatrixXd*);
    void createMatrix();
    inline void addToValue(std::size_t s1, std::size_t g1, std::size_t s2, std::size_t g2, double value)
        { matrixA.coeffRef(idx(s1,g1), idx(s2,g2)) += value; };
    std::function<void(std::size_t, std::size_t, std::size_t, std::size_t, double)> setValue;
};

class Core
{
private:
    MatrixSystem msys;
    Sizing* sz;
    System* sys;
    
    Eigen::MatrixXd independentTerms; // "independent terms" equals b, the right-hand side in Ax=b when solving for x [species_index, grid_location]
    Eigen::MatrixXd gridConcentration, gridConcentrationPrevious; // current concentrations in grid for all species [species_index, grid_location]
    
    // matrix that stores the 2nd order kinetic interactions (needed for Laasonen linearization)
    // secondOrderKinetics[] = {A, B, C, normalized_rate} for the 2nd order reaction A+B<->C
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t, double>> secondOrderKinetics;
    std::vector<std::tuple<std::size_t, std::size_t, double>> firstOrderKinetics;
    std::vector<bool> speciesInRedox;
    std::vector<double> matrixB0DiagonalTerms;
    Eigen::MatrixXd currentContributionMatrix;
    Eigen::VectorXd currentContributionSpeciesFlux, currentContributionRedoxFlux;
public:
    void initialize(Sizing*, System*, double);
    void initVectorsEtc();
    
    void solveSystem(double);
    void invertMatrix();
    void updateIndependentTerms();
    void updateRedoxInMatrix(double);
    void updateKineticsInMatrix(double);
    double calcCurrentFromFlux();

    void addRedoxToMatrix(double);
    void addBICoeffsToMatrix();
    void addKineticsToMatrix();
    void addHalfReactionKineticsToMatrix(Species*, Species*, Species*, Species*, double);
    void addFirstOrderKineticTerm(Species*, Species*, double);
    void addSecondOrderKineticTerm(Species*, Species*, Species*, double);
};

#endif // SIMULATIONCORE_H
