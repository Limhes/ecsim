# code description

(Work in progress!)

**Main files**

`simulation.h/cpp` contains class Simulation, which is the main class you need to instantiate. This class is the interface between the simulation setup (system/species/redox/reaction, electrode, environment, experiment), the simulation algorithm, and the user. Uses *dimensioned* variables (V, s, ...). This class calls the core algorithm at each potential step.

`simulationcore.h/cpp` contains the simulation algorithm. Uses *dimensionless* variables. The Sizing class handles variable normalization and stores the various dimensions used in the algorithm (number of species, reactions, redox steps, derivative/current coefficients, grid points, etc.). The Core class is the simulation algorithm, which sets up the matrix system and then solves it successively at each potential step. The MatrixSystem class is a wrapper around Eigen's SparseMatrix & LU solver which is the real-real core of the algorithm.

**Other files:**

* `system.h/cpp` contains classes System, Species, Redox, Reaction
* `electrodes.h/cpp` contains class Electrode
* `environment.h/cpp` contains class Environment
* `experiment.h/cpp` contains class Experiment
* `coefs_alpha_beta.h/cpp` contains the alpha and beta coefficients and is a direct copy from Molina et al. (10.1016/j.cplett.2015.11.011). Expressions for these coefficients can be found here: 10.1016/j.electacta.2008.08.039
