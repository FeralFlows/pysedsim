""" Python wrapper for the Borg MOEA.

Provides a Python interface for the Borg MOEA.  The Borg MOEA shared library
(typically named libborg.so or borg.dll) must be located in the same directory
as this module.  A simple example of using this module is provided below.

    from borg import *

    borg = Borg(2, 2, 0, lambda x,y : [x**2 + y**2, (x-2)**2 + y**2],
        bounds=[[-50, 50], [-50, 50]],
        epsilons=[0.01, 0.01])

    for solution in borg.solve({'maxEvaluations':10000}):
        solution.display()

This wrapper can also run the master-slave and multi-master implementations
of the Borg MOEA.  

    Configuration.startMPI()
    borg = Borg(...)
    borg.solveMPI(islands=4, maxTime=1)
    Configuration.stopMPI()

Please cite the following paper in any works that use or are derived from this
program.

    Hadka, D. and Reed, P. (2013).  "Borg: An Auto-Adaptive Many-Objective
    Evolutionary Computing Framework."  Evolutionary Computation,
    21(2):231-259.

Copyright 2013-2014 David Hadka
Requires Python 2.5 or later
"""

from ctypes import *
import os
import sys
import time

terminate=False

class Configuration:
    """ Holds configuration options for the Borg MOEA Python wrapper. """

    @staticmethod
    def check():
        """ Checks if the Borg MOEA is initialized and ready to run; otherwise an error is raised. """
        try:
            Configuration.libc
        except:
            raise OSError("The standard C library is not defined, please see Configuration.setStandardCLibrary(<file>)")

        try:
            Configuration.libborg
        except:
            raise OSError("The Borg MOEA C library is not defined, please see Configuration.setBorgLibrary(<file>)")

    @staticmethod
    def initialize():
        """ Initializes the standard C and Borg MOEA libraries. """
        Configuration.setStandardCLibrary()
        Configuration.setBorgLibrary()
        Configuration.seed()
        Configuration.startedMPI = False

    @staticmethod
    def setStandardCLibrary(path=None):
        """ Override the standard C library (libc) used by the Python-to-C interface.

        If the path is not specified, this method will attempt to auto-detect the
        correct location of the standard C library.  If this auto-detection fails,
        this method will return without error.  This allows the module to load
        successfully and requires the user to manually invoke this method before
        using the Borg MOEA.
        """

        if path:
            Configuration.libc = CDLL(path)
        elif os.name == "posix":
            try:
                Configuration.libc = CDLL("libc.so.6")
            except OSError:
                return
        elif os.name == "nt" and cdll.msvcrt:
            Configuration.libc = cdll.msvcrt
        else:
            return

        try:
            Configuration.stdout = Configuration.libc.fdopen(sys.stdout.fileno(), "w")
        except AttributeError:
            Configuration.stdout = Configuration.libc._fdopen(sys.stdout.fileno(), "w")

    @staticmethod
    def setBorgLibrary(path=None):
        """ Override the location of the Borg MOEA shared object.

        If the path is not specified, this method attempts to auto-detect the location
        of the Borg MOEA C library.  If auto-detection fails, this method returns
        without error.  This allows the module to load successfully and requires the
        user to manually invoke this method before using the Borg MOEA
        """

        if path:
            try:
                Configuration.libborg = CDLL(path)
                Configuration.libborg.BORG_Copyright
                Configuration.stdcall = False
            except AttributeError:
                # Not using __cdecl, try __stdcall instead
                if os.name == "nt":
                    Configuration.libborg = WinDLL(path)
                    Configuration.stdcall = True
        elif os.name == "posix":
            try:
                Configuration.libborg = CDLL("./libborg.so")
                Configuration.stdcall = False
            except OSError:
                return
        elif os.name == "nt":
            try:
                Configuration.libborg = CDLL("./borg.dll")
                Configuration.libborg.BORG_Copyright
                Configuration.stdcall = False
            except OSError:
                return
            except AttributeError:
                # Not using __cdecl, try __stdcall instead
                try:
                    Configuration.libborg = WinDLL("./borg.dll")
                    Configuration.stdcall = True
                except OSError:
                    return
                    
        # Set result type of functions with non-standard types
        Configuration.libborg.BORG_Solution_get_variable.restype = c_double
        Configuration.libborg.BORG_Solution_get_objective.restype = c_double
        Configuration.libborg.BORG_Solution_get_constraint.restype = c_double
        Configuration.libborg.BORG_Operator_get_probability.restype = c_double

    @staticmethod
    def seed(value=None):
        """ Sets the pseudo-random number generator seed. """
        Configuration.check()

        if value:
            Configuration.libborg.BORG_Random_seed(c_ulong(value))
        else:
            Configuration.libborg.BORG_Random_seed(c_ulong(os.getpid()*long(time.time())))

    @staticmethod
    def enableDebugging():
        """ Enables debugging output from the Borg MOEA. """
        Configuration.check()
        Configuration.libborg.BORG_Debug_on()

    @staticmethod
    def disableDebugging():
        """ Disables debugging output from the Borg MOEA. """
        Configuration.check()
        Configuration.libborg.BORG_Debug_off()

    @staticmethod
    def displayCopyright():
        """ Displays the copyright message for the Borg MOEA. """
        Configuration.check()
        Configuration.libborg.BORG_Copyright(Configuration.stdout)

    @staticmethod
    def startMPI():
        """ Initializes MPI to enable master-slave and multi-master Borg MOEA runs. """
        if Configuration.startedMPI:
            raise RuntimeError("MPI is already started")

        if os.name != "posix":
            raise RuntimeError("MPI is only supported on Linux")

        try:
            Configuration.libborg.BORG_Algorithm_ms_startup
        except AttributeError:
            # The serial Borg MOEA C library is loaded; switch to parallel
            try:
                Configuration.setBorgLibrary("./libborgmm.so")
            except OSError:
                try:
                    Configuration.setBorgLibrary("./libborgms.so")
                except OSError:
                    raise OSError("Unable to locate the parallel Borg MOEA C library")

        # The following line is needed to load the MPI library correctly
        CDLL("libmpi.so", RTLD_GLOBAL)

        # Pass the command-line arguments to MPI_Init
        argc = c_int(len(sys.argv))
        CHARPP = c_char_p * len(sys.argv)
        argv = CHARPP()

        for i in range(len(sys.argv)):
            argv[i] = sys.argv[i]

        Configuration.libborg.BORG_Algorithm_ms_startup(
            cast(addressof(argc), POINTER(c_int)),
            cast(addressof(argv), POINTER(CHARPP)))

        Configuration.startedMPI = True

    @staticmethod
    def stopMPI():
        """ Shuts down MPI; the master-slave and multi-master Borg MOEA can no longer be used. """
        if not Configuration.startedMPI:
            raise RuntimeError("MPI is not started")

        Configuration.libborg.BORG_Algorithm_ms_shutdown()
        Configuration.startedMPI = False

class RestartMode:
    """ Controls the mutation rate during restarts.

    DEFAULT  - The mutation rate is fixed at 1/numberOfVariables
    RANDOM   - The mutation rate is fixed at 100%
    RAMPED   - The mutation rates are uniformly sampled between 1/numberOfVariables to 100%
    ADAPTIVE - The mutation rate adapts based on success of previous restarts
    INVERTED - Similar to ADAPTIVE, except the rate is inverted
    """

    DEFAULT  = 0
    RANDOM   = 1
    RAMPED   = 2
    ADAPTIVE = 3
    INVERTED = 4

class ProbabilityMode:
    """ Controls how operator probabilities are adapted.

    DEFAULT  - Operator probabilities based on archive membership
    RECENCY  - Operator probabilities based on recency (tracks recent additions to archive)
    BOTH     - Operator probabilities based on archive membership and recency
    ADAPTIVE - Favors archive membership, but uses recency if insufficient archive size
    """

    DEFAULT  = 0
    RECENCY  = 1
    BOTH     = 2
    ADAPTIVE = 3

class InitializationMode:
    """ Controls how initial populations in the multi-master Borg MOEA are initialized.

    UNIFORM      - Each master starts with a uniformly distributed population
    LATIN        - Each master starts with a Latin hypercube sampled population
    GLOBAL_LATIN - A global Latin hypercube sampled population is generated, partitioned,
               and distributed to the master nodes
    """

    UNIFORM      = 0
    LATIN        = 1
    GLOBAL_LATIN = 2

class Direction:
    """ The optimization direction of an objective (minimized or maximized).

    MINIMIZE - The objective is minimized towards negative infinity
    MAXIMIZE - The objective is maximized towards positive infinity
    """

    MINIMIZE = 0
    MAXIMIZE = 1

class Borg:
    """ Solves an optimization problem using the Borg MOEA. """

    def __init__(self, numberOfVariables, numberOfObjectives, numberOfConstraints, function, epsilons=None,
                 bounds=None, directions=None, add_pysedsim_inputs=None):
        """ Creates a new instance of the Borg MOEA.

        numberOfVariables   - The number of decision variables in the optimization problem
        numberOfObjectives  - The number of objectives in the optimization problem
        numberOfConstraints - The number of constraints in the optimization problem
        function            - The function defining the optimization problem
        epsilons            - The epsilon values for each objective
        bounds              - The lower and upper bounds for each decision variable
        directions          - The optimization direction (MINIMIZE or MAXIMIZE) for each objective
        """

        # Ensure the underlying library is available
        Configuration.check()

        # Validate input arguments
        if numberOfVariables < 1:
            raise ValueError("Requires at least one decision variable")

        if numberOfObjectives < 1:
            raise ValueError("Requires at least one objective")

        if numberOfConstraints < 0:
            raise ValueError("Number of constraints can not be negative")

        # Construct Borg object
        self.numberOfVariables = numberOfVariables
        self.numberOfObjectives = numberOfObjectives
        self.numberOfConstraints = numberOfConstraints
        self.directions = directions
        if add_pysedsim_inputs is None:
            self.function = _functionWrapper(function, numberOfVariables, numberOfObjectives, numberOfConstraints, directions)
        else:
            # More PySedSim inputs are required
            self.function = _functionWrapper(function, numberOfVariables, numberOfObjectives, numberOfConstraints,
                                             directions, addl_inputs=add_pysedsim_inputs)

        if Configuration.stdcall:
            self.CMPFUNC = WINFUNCTYPE(c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double))
        else:
            self.CMPFUNC = CFUNCTYPE(c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double))

        self.callback = self.CMPFUNC(self.function)
        self.reference = c_void_p(Configuration.libborg.BORG_Problem_create(numberOfVariables, numberOfObjectives, numberOfConstraints, self.callback))

        if bounds:
            self.setBounds(*bounds)
        else:
            self.setBounds(*[[0, 1]]*numberOfVariables)

        if epsilons:
            self.setEpsilons(*epsilons)
        else:
            self.epsilonsAssigned = False

    def __del__(self):
        """ Deletes the underlying C objects. """
        try:
            Configuration.libborg.BORG_Problem_destroy(self.reference)
        except AttributeError:
            pass

    def setBounds(self, *args):
        """ Sets the decision variable lower and upper bounds.

        The arguments to this function must be 2-ary lists defining the
        lower and upper bounds.  The number of lists must equal the
        number of decision variables.  For example:
            setBounds([0, 1], [-10, 10], [-1, 1])
        If each decision variable has the same bounds, this can be
        written compactly:
            setBounds(*[[0, 1]]*3)
        """

        if len(args) != self.numberOfVariables:
            raise ValueError("Incorrect number of bounds specified")

        for i in range(self.numberOfVariables):
            self._setBounds(i, args[i][0], args[i][1])

    def setEpsilons(self, *args):
        """ Sets the epsilons for the objective values.

        The epsilons control the granularity / resolution of the Pareto
        optimal set.  Small epsilons typically result in larger Pareto
        optimal sets, but can reduce runtime performance.  Specify one
        argument for each objective.  For example:
            setEpsilons(0.01, 0.5)
        If all epsilons are the same, this can be written more compactly:
            setEpsilons(*[0.01]*2)
        """

        if len(args) != self.numberOfObjectives:
            raise ValueError("Incorrect number of epsilons specified")

        for i in range(self.numberOfObjectives):
            self._setEpsilon(i, args[i])

        self.epsilonsAssigned = True

    def _setEpsilon(self, index, value):
        """ Sets the epsilon value at the given index. """
        Configuration.libborg.BORG_Problem_set_epsilon(self.reference, index, c_double(value))

    def _setBounds(self, index, lowerBound, upperBound):
        """ Sets the lower and upper decision variable bounds at the given index. """
        Configuration.libborg.BORG_Problem_set_bounds(self.reference, index, c_double(lowerBound), c_double(upperBound))

    def solveMPI(self, islands=1, maxTime=None, maxEvaluations=None, initialization=None, runtime=None,
            allEvaluations=None, frequency=None):
        """ Runs the master-slave or multi-master Borg MOEA using MPI.

        islands        - The number of islands
        maxTime        - The maximum wallclock time to run, in hours
        maxEvaluations - The maximum NFE per island (total NFE is islands*maxEvaluations)
        initialization - Controls how the initial populations are generated
        runtime        - Filename pattern for saving runtime dynamics (the filename should include
                 one %d which gets replaced by the island index)
        allEvaluations - Filename pattern for saving all evaluations (the filename should include
                         one %d which gets replaced by the island index).  Since this can quickly
                         generate large files, use this option with caution.
        frequency      - Frequency of runtime output as integer (e.g., 500). Default is 100.

        Note: All nodes must invoke solveMPI.  However, only one node will return the discovered
        Pareto optimal solutions.  The rest will return None.
        """

        if not self.epsilonsAssigned:
            raise RuntimeError("Epsilons must be assigned")

        if not Configuration.startedMPI:
            raise RuntimeError("MPI is not started; call Configuration.startMPI() first")

        if not maxTime and not maxEvaluations:
            raise ValueError("Must specify maxEvaluations or maxTime (or both)")

        if islands > 1:
            try:
                Configuration.libborg.BORG_Algorithm_ms_islands(c_int(islands))
            except AttributeError:
                raise RuntimeError("The loaded Borg MOEA C library does not support multi-master")

        if maxTime:
            Configuration.libborg.BORG_Algorithm_ms_max_time(c_double(maxTime))

        if maxEvaluations:
            Configuration.libborg.BORG_Algorithm_ms_max_evaluations(c_int(maxEvaluations))

        if initialization and islands > 1:
            Configuration.libborg.BORG_Algorithm_ms_initialization(c_int(initialization));

        if runtime:
            Configuration.libborg.BORG_Algorithm_output_runtime(c_char_p(runtime));

        if frequency:
            Configuration.libborg.BORG_Algorithm_output_frequency(c_char_p(frequency));

        if allEvaluations:
            Configuration.libborg.BORG_Algorithm_output_evaluations(c_char_p(allEvaluations));

        result = Configuration.libborg.BORG_Algorithm_ms_run(self.reference)

        return Result(result, self) if result else None

    def solve(self, settings={}):
        """ Runs the Borg MOEA to solve the defined optimization problem, returning the
        discovered Pareto optimal set.

        settings - Dictionary of parameters for the Borg MOEA.  The key should match one
               of the parameters defined by the C Borg API.  Default parameter values
               are used for any undefined parameters.
        """

        if not self.epsilonsAssigned:
            raise RuntimeError("Epsilons must be set")

        maxEvaluations = settings.get("maxEvaluations", 10000)
        start = time.clock()

        pm = Configuration.libborg.BORG_Operator_create("PM", 1, 1, 2, Configuration.libborg.BORG_Operator_PM)
        Configuration.libborg.BORG_Operator_set_parameter(pm, 0, c_double(settings.get("pm.rate", 1.0 / self.numberOfVariables)))
        Configuration.libborg.BORG_Operator_set_parameter(pm, 1, c_double(settings.get("pm.distributionIndex", 20.0)))
        
        sbx = Configuration.libborg.BORG_Operator_create("SBX", 2, 2, 2, Configuration.libborg.BORG_Operator_SBX)
        Configuration.libborg.BORG_Operator_set_parameter(sbx, 0, c_double(settings.get("sbx.rate", 1.0)))
        Configuration.libborg.BORG_Operator_set_parameter(sbx, 1, c_double(settings.get("sbx.distributionIndex", 15.0)))
        Configuration.libborg.BORG_Operator_set_mutation(sbx, pm)

        de = Configuration.libborg.BORG_Operator_create("DE", 4, 1, 2, Configuration.libborg.BORG_Operator_DE)
        Configuration.libborg.BORG_Operator_set_parameter(de, 0, c_double(settings.get("de.crossoverRate", 0.1)))
        Configuration.libborg.BORG_Operator_set_parameter(de, 1, c_double(settings.get("de.stepSize", 0.5)))
        Configuration.libborg.BORG_Operator_set_mutation(de, pm)

        um = Configuration.libborg.BORG_Operator_create("UM", 1, 1, 1, Configuration.libborg.BORG_Operator_UM)
        Configuration.libborg.BORG_Operator_set_parameter(um, 0, c_double(settings.get("um.rate", 1.0 / self.numberOfVariables)))

        spx = Configuration.libborg.BORG_Operator_create("SPX", c_int(settings.get("spx.parents", 10)), c_int(settings.get("spx.offspring", 2)), 1, Configuration.libborg.BORG_Operator_SPX)
        Configuration.libborg.BORG_Operator_set_parameter(spx, 0, c_double(settings.get("spx.epsilon", 3.0)))

        pcx = Configuration.libborg.BORG_Operator_create("PCX", c_int(settings.get("pcx.parents", 10)), c_int(settings.get("pcx.offspring", 2)), 2, Configuration.libborg.BORG_Operator_PCX)
        Configuration.libborg.BORG_Operator_set_parameter(pcx, 0, c_double(settings.get("pcx.eta", 0.1)))
        Configuration.libborg.BORG_Operator_set_parameter(pcx, 1, c_double(settings.get("pcx.zeta", 0.1)))

        undx = Configuration.libborg.BORG_Operator_create("UNDX", c_int(settings.get("undx.parents", 10)), c_int(settings.get("undx.offspring", 2)), 2, Configuration.libborg.BORG_Operator_UNDX)
        Configuration.libborg.BORG_Operator_set_parameter(undx, 0, c_double(settings.get("undx.zeta", 0.5)))
        Configuration.libborg.BORG_Operator_set_parameter(undx, 1, c_double(settings.get("undx.eta", 0.35)))

        algorithm = Configuration.libborg.BORG_Algorithm_create(self.reference, 6)
        Configuration.libborg.BORG_Algorithm_set_operator(algorithm, 0, sbx)
        Configuration.libborg.BORG_Algorithm_set_operator(algorithm, 1, de)
        Configuration.libborg.BORG_Algorithm_set_operator(algorithm, 2, pcx)
        Configuration.libborg.BORG_Algorithm_set_operator(algorithm, 3, spx)
        Configuration.libborg.BORG_Algorithm_set_operator(algorithm, 4, undx)
        Configuration.libborg.BORG_Algorithm_set_operator(algorithm, 5, um)

        Configuration.libborg.BORG_Algorithm_set_initial_population_size(algorithm, c_int(settings.get("initialPopulationSize", 100)))
        Configuration.libborg.BORG_Algorithm_set_minimum_population_size(algorithm, c_int(settings.get("minimumPopulationSize", 100)))
        Configuration.libborg.BORG_Algorithm_set_maximum_population_size(algorithm, c_int(settings.get("maximumPopulationSize", 10000)))
        Configuration.libborg.BORG_Algorithm_set_population_ratio(algorithm, c_double(1.0 / settings.get("injectionRate", 0.25)))
        Configuration.libborg.BORG_Algorithm_set_selection_ratio(algorithm, c_double(settings.get("selectionRatio", 0.02)))
        Configuration.libborg.BORG_Algorithm_set_restart_mode(algorithm, c_int(settings.get("restartMode", RestartMode.DEFAULT)))
        Configuration.libborg.BORG_Algorithm_set_max_mutation_index(algorithm, c_int(settings.get("maxMutationIndex", 10)))
        Configuration.libborg.BORG_Algorithm_set_probability_mode(algorithm, c_int(settings.get("probabilityMode", ProbabilityMode.DEFAULT)))

        runtimeformat = settings.get('runtimeformat', 'optimizedv')
        fp = None
        if "frequency" in settings:
            statistics = []
            lastSnapshot = 0
            frequency = settings.get("frequency")    
            if "runtimefile" in settings:
                fp = open(settings['runtimefile'], 'w')
                if runtimeformat == 'optimizedv':
                    fp.write("//")
                    dynamics_header = [
                            "NFE", "ElapsedTime", 
                            "SBX", "DE", "PCX", "SPX", "UNDX", "UM",
                            "Improvements", "Restarts", 
                            "PopulationSize", "ArchiveSize"]
                    if settings.get("restartMode", None) == RestartMode.ADAPTIVE:
                        dynamics_header.append("MutationIndex")
                    fp.write(",".join(dynamics_header))
                    fp.write("\n")
                    header = ["NFE"] \
                           + ["dv{0}".format(i) for i in range(self.numberOfVariables)] \
                           + ["obj{0}".format(i) for i in range(self.numberOfObjectives)] \
                           + ["con{0}".format(i) for i in range(self.numberOfConstraints)]
                    fp.write(",".join(header))
                    fp.write("\n")
                    fp.flush()
            else:
                fp = None
        else:
            statistics = None

        data_header_written=False
        while Configuration.libborg.BORG_Algorithm_get_nfe(algorithm) < maxEvaluations:
            Configuration.libborg.BORG_Algorithm_step(algorithm)
            if terminate is True:
                    break
            currentEvaluations = Configuration.libborg.BORG_Algorithm_get_nfe(algorithm)

            if statistics is not None and currentEvaluations-lastSnapshot >= frequency:
                entry = {}
                entry["NFE"] = currentEvaluations
                entry["ElapsedTime"] = time.clock() - start
                entry["SBX"] = Configuration.libborg.BORG_Operator_get_probability(sbx)
                entry["DE"] = Configuration.libborg.BORG_Operator_get_probability(de)
                entry["PCX"] = Configuration.libborg.BORG_Operator_get_probability(pcx)
                entry["SPX"] = Configuration.libborg.BORG_Operator_get_probability(spx)
                entry["UNDX"] = Configuration.libborg.BORG_Operator_get_probability(undx)
                entry["UM"] = Configuration.libborg.BORG_Operator_get_probability(um)
                entry["Improvements"] = Configuration.libborg.BORG_Algorithm_get_number_improvements(algorithm)
                entry["Restarts"] = Configuration.libborg.BORG_Algorithm_get_number_restarts(algorithm)
                entry["PopulationSize"] = Configuration.libborg.BORG_Algorithm_get_population_size(algorithm)
                entry["ArchiveSize"] = Configuration.libborg.BORG_Algorithm_get_archive_size(algorithm)

                if settings.get("restartMode", RestartMode.DEFAULT) == RestartMode.ADAPTIVE:
                    entry["MutationIndex"] = Configuration.libborg.BORG_Algorithm_get_mutation_index(algorithm)
                if fp is None:
                    statistics.append(entry)
                else:
                    archive = Result(Configuration.libborg.BORG_Algorithm_get_result(algorithm), self, statistics)
                    if runtimeformat == 'optimizedv':
                        row = ["{0}".format(entry[dynamic]) for dynamic in dynamics_header]
                        fp.write("//")
                        fp.write(",".join(row))
                        fp.write("\n")
                        delimiter = ','
                    elif runtimeformat == 'borg':
                        metrics = [ ("NFE", 'd'),
                            ("ElapsedTime", '.17g'),
                            ("SBX", '.17g'),
                            ("DE", '.17g'),
                            ("PCX", '.17g'),
                            ("SPX", '.17g'),
                            ("UNDX", '.17g'),
                            ("UM", '.17g'),
                            ("Improvements", 'd'),
                            ("Restarts", 'd'),
                            ("PopulationSize", 'd'),
                            ("ArchiveSize", 'd')
                        ]
                        for metric,fmt in metrics:
                            fp.write("//{0}={1}\n".format(metric,"".join(["{0:",fmt,"}"])).format(entry[metric]))
                        if 'MutationIndex' in entry:
                            fp.write("//MutationIndex={0:d}\n".format(entry['MutationIndex']))
                        if "data_header" in settings and data_header_written is False:
                            data_header_written=True
                            data_header = ["_".join(x.split(" ")) for x in settings['data_header']]
                            data_header.insert(0, "NFE")
                            fp.write(" ".join(data_header))
                            fp.write("\n")
                        delimiter = " "

                    for solution in archive:
                        report = [entry["NFE"]]
                        report.extend(solution.getVariables())
                        report.extend(solution.getObjectives())
                        report.extend(solution.getConstraints())
                        fp.write(delimiter.join("{0}".format(v) for v in report))
                        fp.write("\n")
                    fp.flush()

                lastSnapshot = currentEvaluations

        result = Configuration.libborg.BORG_Algorithm_get_result(algorithm)
        if "runtimefile" in settings:
            fp.close()

        Configuration.libborg.BORG_Operator_destroy(sbx)
        Configuration.libborg.BORG_Operator_destroy(de)
        Configuration.libborg.BORG_Operator_destroy(pm)
        Configuration.libborg.BORG_Operator_destroy(um)
        Configuration.libborg.BORG_Operator_destroy(spx)
        Configuration.libborg.BORG_Operator_destroy(pcx)
        Configuration.libborg.BORG_Operator_destroy(undx)
        Configuration.libborg.BORG_Algorithm_destroy(algorithm)

        return Result(result, self, statistics)

class Solution:
    """ A solution to the optimization problem. """

    def __init__(self, reference, problem):
        """ Creates a solution given a reference to the underlying C object. """
        self.reference = reference
        self.problem = problem

    # There is no __del__ since the underlying C solutions are deleted when the associated
    # result object is deleted 

    def getVariables(self):
        """ Returns the decision variable values for this solution. """
        return [self._getVariable(i) for i in range(self.problem.numberOfVariables)]

    def getObjectives(self):
        """ Returns the objective values for this solution. """
        return [self._getObjective(i) for i in range(self.problem.numberOfObjectives)]

    def getConstraints(self):
        """ Returns the constraint values for this solution. """
        return [self._getConstraint(i) for i in range(self.problem.numberOfConstraints)]

    def _getVariable(self, index):
        """ Returns the decision variable at the given index. """
        return Configuration.libborg.BORG_Solution_get_variable(self.reference, index)

    def _getObjective(self, index):
        """ Returns the objective value at the given index. """
        value = Configuration.libborg.BORG_Solution_get_objective(self.reference, index)
    
        if self.problem.directions and self.problem.directions[index]:
            return -value
        else:
            return value

    def _getConstraint(self, index):
        """ Returns the constraint value at the given index. """
        return Configuration.libborg.BORG_Solution_get_constraint(self.reference, index)

    def display(self, out=sys.stdout, separator=" "):
        """ Prints the decision variables, objectives, and constraints to standard output. """
        print >> out, separator.join(map(str, self.getVariables() + self.getObjectives() + self.getConstraints()))

    def violatesConstraints(self):
        """ Returns True if this solution violates one or more constraints; False otherwise. """
        return Configuration.libborg.BORG_Solution_violates_constraints(self.reference) != 0

class Result:
    """ A Pareto optimal set (the output of the Borg MOEA). """

    def __init__(self, reference, problem, statistics=None):
        """ Creates a new Pareto optimal set given a reference to the underlying C object. """
        self.reference = reference
        self.problem = problem
        self.statistics = statistics

    def __del__(self):
        """ Deletes the underlying C objects. """
        Configuration.libborg.BORG_Archive_destroy(self.reference)

    def __iter__(self):
        """ Returns an iterator over the Pareto optimal solutions. """
        return ResultIterator(self)

    def display(self, out=sys.stdout, separator=" "):
        """ Print the Pareto optimal solutions to standard output. """
        for solution in self:
            solution.display(out, separator)

    def size(self):
        """ Returns the size of the Pareto optimal set. """
        return Configuration.libborg.BORG_Archive_get_size(self.reference)

    def get(self, index):
        """ Returns the Pareto optimal solution at the given index. """
        return Solution(Configuration.libborg.BORG_Archive_get(self.reference, index), self.problem)

class ResultIterator:
    """ Iterates over the solutions in a Pareto optimal set. """

    def __init__(self, result):
        """ Creates an iterator over the given Pareto optimal set. """
        self.result = result
        self.index = -1

    def next(self):
        """ Returns the next Pareto optimal solution in the set. """
        self.index = self.index + 1

        if self.index >= self.result.size():
            raise StopIteration
        else:
            return self.result.get(self.index)

def _functionWrapper(function, numberOfVariables, numberOfObjectives, numberOfConstraints, directions=None,
                     addl_inputs=None):
    """ Wraps a Python evaluation function and converts it to the function signature
    required by the C API.

    function - The Python evaluation function of the form (o, c) = f(v)
    numberOfVariables - The number of decision variables
    numberOfObjectives - The number of objectives
    numberOfConstraints - The number of constraints
    directions - The array of optimization directions
    """

    def innerFunction(v,o,c):
        """ The function that gets passed to the C API.

        v - The array of decision variables (input)
        o - The array of objectives (output)
        c - The array of constraint values (output)
        """
        global terminate
        try:
            if addl_inputs is None:
                result = function(*[v[i] for i in range(numberOfVariables)])
            else:
                result = function([v[i] for i in range(numberOfVariables)], addl_inputs)

            objectives = None
            constraints = None

            if isinstance(result, tuple):
                if len(result) > 0:
                    objectives = result[0]
                if len(result) > 1:
                    constraints = result[1]
            elif isinstance(result, list):
                objectives = result
            else:
                objectives = [result]

            if objectives:
                if len(objectives) != numberOfObjectives:
                    raise ValueError("Incorrect number of objectives returned by function")
                for i in range(len(objectives)):
                    if directions and directions[i]:
                        o[i] = -objectives[i]
                    else:
                        o[i] = objectives[i]
            elif numberOfObjectives > 0:
                raise ValueError("No objectives returned by function")

            if constraints:
                if len(constraints) != numberOfConstraints:
                    raise ValueError("Incorrect number of constraints returned by function")
                for i in range(len(constraints)):
                    c[i] = constraints[i]
            elif numberOfConstraints > 0:
                raise ValueError("No constraints returned by function")

            return 0
        except KeyboardInterrupt:
            terminate=True
            return 1
    return innerFunction

class Constraint:
    """ Helper functions for defining constraints.

    These functions ensure several conditions hold.  First, if the
    constraint is satisfied, the value is 0.  If the constraint is
    violated, then the value is non-zero and will scale linearly
    with the degree of violation.
    """

    precision = 0.1

    @staticmethod
    def greaterThan(x, y, epsilon=0.0):
        """ Defines the constraint x > y. """
        return 0.0 if x > y-epsilon else y-x+Constraint.precision

    @staticmethod
    def lessThan(x, y, epsilon=0.0):
        """ Defines the constraint x < y. """
        return 0.0 if x < y+epsilon else x-y+Constraint.precision

    @staticmethod
    def greaterThanOrEqual(x, y, epsilon=0.0):
        """ Defines the constraint x >= y. """
        return 0.0 if x >= y-epsilon else y-x+Constraint.precision

    @staticmethod
    def lessThanOrEqual(x, y, epsilon=0.0):
        """ Defines the constraint x <= y. """
        return 0.0 if x <= y+epsilon else x-y+Constraint.precision

    @staticmethod
    def equal(x, y, epsilon=0.0):
        """ Defines the constraint x == y. """
        return 0.0 if abs(y-x) < epsilon else abs(y-x)+Constraint.precision

    @staticmethod
    def zero(x, epsilon=0.0):
        """ Defines the constraint x == 0. """
        return Constraint.equal(x, 0.0, epsilon)

    @staticmethod
    def nonNegative(x, epsilon=0.0):
        """ Defines the constraint x >= 0. """
        return Constraint.greaterThanOrEqual(x, 0.0, epsilon)

    @staticmethod
    def positive(x, epsilon=0.0):
        """ Defines the constraint x > 0. """
        return Constraint.greaterThan(x, 0.0, epsilon)

    @staticmethod
    def negative(x, epsilon=0.0):
        """ Defines the constraint x < 0. """
        return Constraint.lessThan(x, 0.0, epsilon)

    @staticmethod
    def all(*args):
        """ Requires all conditions to be satisfied. """
        return sum(args)

    @staticmethod
    def any(*args):
        """ Requres at least one condition to be satisfied. """
        return 0.0 if 0.0 in args else sum(args)

Configuration.initialize()
