#ifndef SIMULATIONCORE_H
#define SIMULATIONCORE_H

#include <vector>
#include <cmath>

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/Core>

using namespace std;

using Triplet = Eigen::Triplet<double>; //triplets (row, column, value) to fill sparse matrix
using Solver  = Eigen::SparseLU<Eigen::SparseMatrix<double>>;  //Turns out that SparseLU is actually fast enough, but only when running as Release :)

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
    double f; // F/RT
    
    // dimensionless parameters:
    double deltaTheta = 0.2; // potential step
    double deltaX; // smallest distance
    double paramGamma = 1.67; // grid expansion factor
    double paramR0, electrodeGeometryFactor; // electrode parameters
    double paramLambda; // kinetic parameter of the system
    double minF = 2.2, maxF = 6.2, minLogRate = 3.0, maxLogRate = 7.0; // parameters that relate to homogeneous reaction rates
    
    size_t numSpecies, numRedox, numReactions;
    size_t numGridPoints, numCurrentPoints = 5, numDerivPoints = 6; // number of points in grid, number of points for calculation of current, number of points for diffusion

    double currentFromFlux; // to calculate current from species' flux
    
    double initialize(System*, Electrode*, Environment*, Experiment*, ostream &);
};

class MatrixSystem
{
private:
    Sizing* sz;
    
    vector<Triplet> tripletContainerMatrixA; // this container holds (row, col, value) triplets before matrixA is assembled
    Eigen::SparseMatrix<double> matrixA; // matrix to be inverted: Ax=b corresponds to matrixA * vecx = vecb
    bool matrixPatternAnalyzed;
    Solver sparseMatrixSolver; // solver that inverts the sparse matrixA using LU factorization
    Eigen::VectorXd vecb, vecx; // vectors, defininition see above
public:
    void initialize(Sizing*);
    
    void solve(Eigen::MatrixXd*, Eigen::MatrixXd*);
    void addToMatrix(size_t, size_t, size_t, size_t, double); // add to tripletcontainer
    void createMatrix(); // create matrix from tripletcontainer
    void changeValue(size_t, size_t, size_t, size_t, double); // change directly in matrix:   matrix.coeffref(x,y) = val
    void addToValue(size_t, size_t, size_t, size_t, double); // sum value directly in matrix: matrix.coeffref(x,y) += val
};

class Core
{
private:
    MatrixSystem msys;
    Sizing* sz;
    System* sys;
    
    vector<double> gridCoordinate, paramGammai, paramGamma2i; // grid coordinates and expansion factor (and its ^2), pre-calculated for all points in grid
    vector<double> coeffAlpha, coeffBeta, coeffBeta0; // coefficients for diffusion (alpha and beta) and current (beta0)
    Eigen::MatrixXd independentTerms; // "independent terms" equals b, the right-hand side in Ax=b when solving for x [species_index, grid_location]
    Eigen::MatrixXd gridConcentration, gridConcentrationPrevious; // current concentrations in grid for all species [species_index, grid_location]
    
    // matrix that stores the 2nd order kinetic interactions (needed for Laasonen linearization)
    // secondOrderKinetics[] = {A, B, C, normalized_rate} for the 2nd order reaction A+B<->C
    vector<tuple<size_t, size_t, size_t, double>> secondOrderKinetics;
    vector<tuple<size_t, size_t, double>> firstOrderKinetics;
    vector<bool> speciesInRedox;
    vector<double> matrixB0DiagonalTerms;
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

    double coeffMatrixN2(size_t, int, double);
};

#endif // SIMULATIONCORE_H
