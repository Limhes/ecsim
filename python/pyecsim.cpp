#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
namespace py = pybind11;
#include <stdlib.h>
#include "../source/simulation.h"

/** 
 * The NullStream class acts as a "silent" version of std::cout, akin /dev/null
 * source: https://stackoverflow.com/questions/8243743/is-there-a-null-stdostream-implementation-in-c-or-libraries
 */
class NullStream : public std::ostream {
    class NullBuffer : public std::streambuf {
    public:
        int overflow(int c) { return c; }
    } m_nb;
public:
    NullStream() : std::ostream(&m_nb) {}
};

/** 
 * The Python module to C++ binding definitions
 */
PYBIND11_MODULE(pyecsim, m) {
    py::options options;
    options.disable_function_signatures();
    
    /** 
    * Species
    */
    py::class_<Species>(m, "Species", R"mydelimiter(
    __init__(name:str, conc:float, diff:float) -> pyecsim.Species
    Create chemical species with name, concentration and diffusion coefficient

    :param name: name of the species (can be anything)
    :type name: str
    :param conc: initial concentration [mol/m^3 or mM]
    :type conc: float
    :param diff: diffusion coefficient [m^2/s]
    :type diff: float
    :return: newly generated species
    :rtype: pyecsim.Species
)mydelimiter")
        .def(py::init<std::string, double, double>())
        .def_readonly("normConc", &Species::normalizedConcentration, "Returns the concentration normalized to all species' concentrations")
        .def_readonly("normEquilConc", &Species::normalizedAndEquilibratedConcentration, "Returns the normalized concentration after equilibration")
        .def_readonly("normDiffCoeff", &Species::normalizedDiffusionConstant, "Returns the diffusion coefficient normalized to all species' coefficients");
    
    /** 
    * Reaction
    */
    py::class_<Reaction>(m, "Reaction", R"mydelimiter(
    __init__(A:pyecsim.Species, B:pyecsim.Species, C:pyecsim.Species, D:pyecsim.Species, k_f:float, k_b:float) -> pyecsim.Reaction
    Create Reaction object for reaction A + B <-> C + D
    
    :param A: left-hand side species 1
    :type A: pyecsim.Species
    :param B: left-hand side species 2 (bimolecular forward reaction) or None (unimolecular forward reaction))
    :type B: pyecsim.Species
    :param C: right-hand side species 1
    :type C: pyecsim.Species
    :param D: right-hand side species 2 (bimolecular backward reaction) or None (unimolecular backward reaction)
    :type D: pyecsim.Species
    :param k_f: forward rate constant with unit [1/s] if type(B)==None or [L/mol/s] if type(B)==Species
    :type k_f: float
    :param k_b: backward rate constant with unit [1/s] if type(D)==None or [L/mol/s] if type(D)==Species
    :type k_b: float
    :return: newly generated homogeneous reaction
    :rtype: pyecsim.Reaction
)mydelimiter")
        .def(py::init<Species*, Species*, Species*, Species*, double, double>())
        .def_readwrite("rateConstantForward", &Reaction::rateConstantForward, "The forward rate constant [1/s or L/mol/s]")
        .def_readwrite("rateConstantBackward", &Reaction::rateConstantBackward, "The backward rate constant [1/s or L/mol/s]")
        .def("enable", [](Reaction &rxn) { rxn.enabled = true; return rxn; }, "Enables the reaction and returns the modified Reaction object")
        .def("disable", [](Reaction &rxn) { rxn.enabled = false; return rxn; }, "Disables the reaction and returns the modified Reaction object");
    
    /** 
    * Redox
    */
    py::class_<Redox>(m, "Redox", R"mydelimiter(
    __init__(ox:pyecsim.Species, red:pyecsim.Species, n_e:int, E:float, k_e:float, alpha:float) -> pyecsim.Redox
    Create Redox object for reaction Ox + e <-> Red with Butler-Volmer kinetics

    :param ox: the oxidized form of the species
    :type ox: pyecsim.Species
    :param red: the reduced form of the species
    :type red: pyecsim.Species
    :param n_e: number of electrons in the electron transfer step [-] (n_e >= 1)
    :type n_e: int
    :param E: standard potential of the electron transfer step [V]
    :type E: float
    :param k_e: heterogeneous rate of the electron transfer step [m/s] 
    :type k_e: float
    :param alpha: symmetry parameter of the electron transfer step [-] (0.0 <= alpha <= 1.0)
    :type alpha: float
    :return: newly generated electron transfer step
    :rtype: pyecsim.Redox
)mydelimiter")
        .def(py::init([](Species* ox, Species* red, int n_e, double E, double k_e, double alpha) { return new Redox(ox, red, n_e, E, k_e, alpha, false); }))
        .def_readwrite("n", &Redox::numberElectrons, "Number of electrons [-]")
        .def_readwrite("E", &Redox::standardPotential, "Standard potential [V]")
        .def_readwrite("k_e", &Redox::rateConstantHetero, "Heterogeneous rate constant [m/s]")
        .def_readwrite("alpha", &Redox::alpha, "Symmetry parameter [-] (0.0 <= alpha <= 1.0)")
        .def("enable", [](Redox &rdx) { rdx.enabled = true; return rdx; }, "Enables the electron transfer step and returns the modified Redox object")
        .def("disable", [](Redox &rdx) { rdx.enabled = false; return rdx; }, "Disables the electron transfer step and returns the modified Redox object");

    /** 
    * System
    */
    py::class_<System>(m, "System", R"mydelimiter(
    __init__() -> pyecsim.System
    Create System object containing all species, reactions and redox steps.
    
    Before the simulation starts, the system is brought into a state of chemical equilibrium.

    :return: newly generated chemical/electrochemical system
    :rtype: pyecsim.System
)mydelimiter")
        .def(py::init<>())
        .def("addSpecies", [](System &sys, Species* spec)
            { sys.addSpecies(spec); return sys; }, R"mydelimiter(
    addSpecies(spec:pyecsim.Species) -> pyecsim.System
    Add single species to the system
    
    :param spec: species to be added to the system
    :type spec: pyecsim.Species
    :return: the modified system
    :rtype: pyecsim.System
)mydelimiter")
        .def("addSpecies", [](System &sys, std::vector<Species*> spec_list)
            { for (auto spec : spec_list) { sys.addSpecies(spec); } return sys; }, R"mydelimiter(
    addSpecies(spec_list:List[pyecsim.Species]) -> pyecsim.System
    Add multiple species to the system

    :param spec_list: list of species to be added to the system
    :type spec_list: List[pyecsim.Species]
    :return: the modified system
    :rtype: pyecsim.System
)mydelimiter")
        .def("addRedox", [](System &sys, Redox* red)
            { sys.addRedox(red); return sys; }, R"mydelimiter(
    addRedox(red:pyecsim.Redox) -> pyecsim.System
    Add single electron transfer step to the system

    :param red: electron transfer step to be added to the system
    :type red: pyecsim.Redox
    :return: the modified system
    :rtype: pyecsim.System
)mydelimiter")
        .def("addRedox", [](System &sys, std::vector<Redox*> red_list)
            { for (auto red : red_list) { sys.addRedox(red); } return sys; }, R"mydelimiter(
    addRedox(red_list:List[pyecsim.Redox]) -> pyecsim.System
    Add multiple electron transfer steps to the system

    :param red_list: list of electron transfer steps to be added to the system
    :type red_list: List[pyecsim.Redox]
    :return: the modified system
    :rtype: pyecsim.System
)mydelimiter")
        .def("addReaction", [](System &sys, Reaction* rxn)
            { sys.addReaction(rxn); return sys; }, R"mydelimiter(
    addReaction(rxn:pyecsim.Reaction) -> pyecsim.System
    Add single reaction to the system

    :param rxn: reaction to be added to the system
    :type rxn: pyecsim.Reaction
    :return: the modified system
    :rtype: pyecsim.System
)mydelimiter")
        .def("addReaction", [](System &sys, std::vector<Reaction*> rxn_list)
            { for (auto rxn : rxn_list) { sys.addReaction(rxn); } return sys; }, R"mydelimiter(
    addReaction(rxn_list:List[pyecsim.Reaction]) -> pyecsim.System
    Add multiple reactions to the system

    :param rxn_list: list of reactions to be added to the system
    :type rxn_list: List[pyecsim.Reaction]
    :return: the modified system
    :rtype: pyecsim.System
)mydelimiter")
        .def("getSpeciesByName", &System::getSpeciesByName, R"mydelimiter(
    getSpeciesByName(name:str) -> pyecsim.Species
    Get species reference by name

    :param name: name of the species
    :type name: String
    :return: reference to the species
    :rtype: pyecsim.Species
)mydelimiter")
        .def("generateUniqueSpeciesName", &System::generateUniqueSpeciesName, R"mydelimiter(
    generateUniqueSpeciesName() -> str
    Generate unique species name

    :return: a species name that is guaranteed to be not present in current system
    :rtype: str
)mydelimiter")
        .def("isSpeciesPresentWithName", &System::isSpeciesPresentWithName, R"mydelimiter(
    isSpeciesPresentWithName(name:str) -> bool
    Check whether species with this name exists

    :param name: name of the species
    :type name: String
    :return: True if species with name is present in system
    :rtype: bool
)mydelimiter");

    /** 
    * Experiment
    */
    py::class_<Experiment>(m, "Experiment", "Experimental potential/time parameters", R"mydelimiter(
    __init__() -> pyecsim.Experiment
    Create voltammetry (LSV & CV) Experiment object (all values are initialized as 0.0)
    
    A voltammetry experiment is run as the following sequence:
    
    * If conditioning_time > 0: hold potential at conditioning_potential [V] for conditioning_time [s]
    * If equilibration_time > 0: apply no potential for equilibration_time [s]
    * Repeat num_cycles times (and record current in the last cycle only):
        * Set potential to initial_potential [V]
        * For each vertex_potential in vertex_potentials:
            * Scan to vertex_potential [V] at scan_rate [V/s]
            * Hold potential at vertex_potential [V] for vertex_delay [s]
        * Scan to final_potential [V] at scan_rate [V/s]
    
    :return: newly generated experiment
    :rtype: pyecsim.Experiment
    
    .. note:: Note that this implies that if initial_potential is not equal to final_potential and num_cycles > 1, the potential jumps, which can lead to non-physical behaviour.
)mydelimiter")
        .def(py::init<>())
        .def("setConditioning", [](Experiment &exp, double time, double potential)
            { exp.setConditioningTime(time); exp.setConditioningPotential(potential); return exp; }, R"mydelimiter(
    setConditioning(ct:float, cp:float)
    Set electrode conditioning time and potential

    :param ct: conditioning time [s]
    :type ct: float
    :param cp: conditioning potential [V]
    :type cp: float
)mydelimiter")
        .def("setEquilibration", &Experiment::setEquilibration, R"mydelimiter(
    setEquilibration(et:float)
    Set equilibration time

    :param et: equilibration time [s]
    :type et: float
)mydelimiter")
        .def("setScanPotentials", [](Experiment &exp, double ip, std::vector<double> vp, double fp)
            { exp.setInitialPotential(ip); exp.setVertexPotentials(vp); exp.setFinalPotential(fp); return exp; }, R"mydelimiter(
    setScanPotentials(ip:float, vp:List[float], fp:float)
    Set CV scan potentials

    :param ip: initial potential [V]
    :type ip: float
    :param vp: vertex potentials [V], pass empty list for linear sweep voltammetry
    :type vp: List[float]
    :param fp: final potential [V]
    :type fp: float
)mydelimiter")
        .def("setVertexDelay", &Experiment::setVertexDelay, R"mydelimiter(
    setVertexDelay(vd:float)
    Set vertex delay time for vertices in cyclic voltammetry scan 

    :param vd: vertex delay [s]
    :type vd: float
)mydelimiter")
        .def("setScanRate", &Experiment::setScanRate, R"mydelimiter(
    setScanRate(sr:float)
    Set scan rate of voltammetry scan
    
    :param sr: scan rate  [V/s]
    :type sr: float
)mydelimiter")
        .def("setNumCycles", &Experiment::setNumCycles, R"mydelimiter(
    setNumCycles(nc:int)
    Set number of CV scan cycles
    
    Current is only recorded on the last scan
    
    :param nc: number of repetitive voltammetry scans [-] (nc >= 1)
    :type nc: int
)mydelimiter")
        .def("totalTime", &Experiment::totalTime, R"mydelimiter(
    totalTime() -> float
    Get total experiment time
    
    :return: total duration of voltammetry experiment including conditioning and equilibration [s]
    :rtype: float
)mydelimiter");
    
    /** 
    * Environment
    */
    py::class_<Environment>(m, "Environment", R"mydelimiter(
    __init__() -> pyecsim.Environment
    Create Environment object

    :return: newly generated physical environment
    :rtype: pyecsim.Environment
)mydelimiter")
        .def(py::init<>())
        .def("setTemperature", &Environment::setTemperature, R"mydelimiter(
    setTemperature(temp:float)
    Set solution temperature
    
    :param temp: temperature [K]
    :type temp: float
)mydelimiter");
        
    /** 
    * Electrode
    */
    py::class_<Electrode>(m, "Electrode", R"mydelimiter(
    __init__() -> pyecsim.Electrode
    Create Electrode object

    :return: newly generated electrode
    :rtype: pyecsim.Electrode
)mydelimiter")
        .def(py::init<>())
        .def("disk", [](Electrode &el, double radius) { el.setType(0); el.setGeom1(radius); return el; }, R"mydelimiter(
    disk(radius:float) -> pyecsim.Electrode
    Set disk electrode
    
    :param radius: radius [m]
    :type radius: float
    :return: modified electrode
    :rtype: pyecsim.Electrode
)mydelimiter")
        .def("square", [](Electrode &el, double width) { el.setType(1); el.setGeom1(width); return el; }, R"mydelimiter(
    square(width:float) -> pyecsim.Electrode
    Set square electrode
    
    :param width: width [m]
    :type width: float
    :return: modified electrode
    :rtype: pyecsim.Electrode
)mydelimiter")
        .def("rectangle", [](Electrode &el, double width, double height) { el.setType(2); el.setGeom1(width); el.setGeom2(height); return el; }, R"mydelimiter(
    rectangle(width:float, height:float) -> pyecsim.Electrode
    Set rectangular electrode
    
    :param width: width [m]
    :type width: float
    :param height: height [m]
    :type height: float
    :return: modified electrode
    :rtype: pyecsim.Electrode
)mydelimiter")
        .def("cylinder", [](Electrode &el, double radius, double length) { el.setType(3); el.setGeom1(radius); el.setGeom2(length); return el; }, R"mydelimiter(
    cylinder(radius:float, length:float) -> pyecsim.Electrode
    Set cylindrical electrode
    
    :param radius: radius [m]
    :type radius: float
    :param length: length [m]
    :type length: float
    :return: modified electrode
    :rtype: pyecsim.Electrode
)mydelimiter")
        .def("sphere", [](Electrode &el, double radius) { el.setType(4); el.setGeom1(radius); return el; }, R"mydelimiter(
    sphere(radius:float) -> pyecsim.Electrode
    Set spherical electrode
    
    :param radius: radius [m]
    :type radius: float
    :return: modified electrode
    :rtype: pyecsim.Electrode
)mydelimiter")
        .def("hemisphere", [](Electrode &el, double radius) { el.setType(5); el.setGeom1(radius); return el; }, R"mydelimiter(
    hemisphere(radius:float) -> pyecsim.Electrode
    Set hemispherical electrode
    
    :param radius: radius [m]
    :type radius: float
    :return: modified electrode
    :rtype: pyecsim.Electrode
)mydelimiter");
    
    /** 
    * Simulation
    */
    py::class_<Simulation>(m, "Simulation", R"mydelimiter(
    __init__(verbose:bool) -> pyecsim.Simulation
    Create Simulation object

    * The core of the simulation is an implementation of Molina et al. "Brute force (or not so brute) digital simulation in electrochemistry revisited" Chemical Physics Letters 643 (2016) 71–76 (DOI: 10.1016/j.cplett.2015.11.011)
    * Various practicalities were implemented according to: Britz "Digital Simulation in Electrochemistry" and Compton "Understanding Voltammetry - Simulation of Electrode Processes"
    * Second-order homogeneous reactions are implemented as Laasonen linearizations (Britz p. 165+)
    * Species concentrations are equilibrated before the simulation starts
    * Electron transfer is implemented as Butler-Volmer kinetics
    * Current contributions are solved in a generic way by writing Compton p. 117-119 in matrix form
    
    :param verbose: False to suppress all output
    :type verbose: bool
    :return: newly generated simulation
    :rtype: pyecsim.Simulation
)mydelimiter")
        .def(py::init([](bool verbose) { return new Simulation( verbose ? std::cout : static_cast<std::ostream&>(*(new NullStream())) ); })) // custom constructor
        .def_readonly("env", &Simulation::env, "In-built Environment object")
        .def_readonly("el", &Simulation::el, "In-built Electrode object")
        .def_readonly("exper", &Simulation::exper, "In-built Experiment object")
        .def_readonly("sys", &Simulation::sys, "In-built System object")
        .def("setGridSizing", &Simulation::setGridSizing, R"mydelimiter(
    setGridSizing(gamma_e:float, minF:float, maxF:float, minlograte:float, maxlograte:float)
    Set exponential grid parameters
    
    cf. the Molina paper, Britz and Compton books and the settings tab on http://limhes.net/ecsim/
    
    :param gamma: grid expansion factor in x_i = h * (gamma_e^i - 1) / (gamma_e - 1) with gamma_e > 1.0
    :type gamma: Float
    :param minF: see the settings tab on http://limhes.net/ecsim/
    :type minF: Float
    :param maxF: see the settings tab on http://limhes.net/ecsim/
    :type maxF: Float
    :param minlograte: see the settings tab on http://limhes.net/ecsim/
    :type minlograte: Float
    :param maxlograte: see the settings tab on http://limhes.net/ecsim/
    :type maxlograte: Float
    
    .. warning:: For advanced users only
)mydelimiter")
        .def("setPotentialSizing", &Simulation::setPotentialSizing, R"mydelimiter(
    setPotentialSizing(deltatheta:float)
    Set dimensionless potential step
    
    cf. Britz and Compton books and the settings tab on http://limhes.net/ecsim/
    
    :param deltatheta: dimensionless potential step
    :type deltatheta: Float
    
    .. warning:: For advanced users only
)mydelimiter")
        .def("setDifferentialOrders", &Simulation::setDifferentialOrders, R"mydelimiter(
    setDifferentialOrders(num_curr:int, num_diff:int)
    Set differential orders for current and diffusion calculations
    
    cf. Molina paper
    
    :param num_curr: number of points for current differential
    :type num_curr: Integer
    :param num_diff: number of points for diffusion differential
    :type num_diff: Integer
    
    .. warning:: For advanced users only
)mydelimiter")
        .def("run", [](Simulation &sim) { std::vector<std::vector<double>> results{{}, {}}; sim.run(results[1], results[0]); return results; }, R"mydelimiter(
    run()
    Run the simulation
    
    :return: the voltammetry simulation potentials and corresponding currents as [potentials, currents]
    :rtype: List of List of Float
)mydelimiter")
        .def("metrics", [](Simulation &sim) { std::vector<double> metrics{sim.ipa, sim.Epa, sim.ipc, sim.Epc}; return metrics; }, R"mydelimiter(
    metrics() -> List[float]
    Returns metrics recorded during the simulation
    
    :return: maximum and minimum currents recorded and the corresponding potentials [max. current, potential at max. current, min. current, potential at min. current]
    :rtype: List[float]
)mydelimiter");
    
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
