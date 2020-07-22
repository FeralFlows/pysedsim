# PySedSim
An Open Source Reservoir and Sediment Simulation Screening Framework for Identifying and Evaluating Dam Siting, Design, and Operation Alternatives.

![sedsim_worksheet](/images/pysedsim_logo.png)

- [Contact Us](#Contact)
- [Introduction](#Introduction)
- [Getting Started](#InstallGuide)
- [Overview of Model File Structure](#ModelFileStructure)
- [Publications using SedSim](#Pubs)

# <a name="Contact Us"></a>Contact
Thomas B. Wild, twild@umd.edu

Patrick M. Reed, patrick.reed@cornell.edu

Abigail N. Birnbaum, birnbaum.abigail@gmail.com

Daniel P. Loucks, loucks@cornell.edu

# <a name="Introduction"></a>Introduction

PySedSim is an open-source, daily time step stochastic river basin simulation model for water and sediment flows, and hydropower production, in networks of reservoirs and river channels. PySedSim enables water resources systems analysts and planners to explore alternative system configurations of reservoir sites, designs (i.e., dam outlet structures), and operating policies (SDO), and their implications for water flows, sediment transport, reservoir sediment trapping, and hydropower production in any river basin. The model enables simulation of a wide range of reservoir sediment management techniques, including flushing, sluicing, density current venting, bypassing, and dredging. The model performs a daily time-step mass-balance simulation of flow and sediment that is intended to predict in relative terms the spatial and temporal accumulation and depletion of sediment in river reaches and in reservoirs under different reservoir operating and sediment management policies. Thus, the model is expected to be used for estimating sediment transport in river basis including those that have experienced (or will experience) extensive reservoir development. PySedSim is generic model, in that it is comprised of flexible code that will let you model any network of river channels and reservoirs of interest to you. You provide PySedSim the basic information about the system you want to model in the form of a handful of .xlsx/.csv input files. PySedSim then creates and simulates the mass balance of water and sediment in that user-specified specified network, according to various preferences. It also exports and stores the results of your simulation in .csv file(s).

The source code is written in Python. PySedSim includes four core model features: 1) representing alternative reservoir sediment management approaches, 2) representing detailed dam design features (e.g., gates), 3) supporting emerging multi-objective optimization frameworks for discovering reservoir operating policies and their resulting tradeoffs, and 4) facilitating stochastic Monte Carlo simulation for characterizing uncertainty in hydrologic and sediment processes. The model strongly broadens the SDO design analytic capabilities of the sediment simulation SedSim model (Wild and Loucks 2015a; Wild et al., 2019)), while providing backward compatibility with its EXCEL-based interfacing for training of less advanced users. The model was developed at Cornell University, in partnership with the Natural Heritage Institute (NHI), as well as at the University of Maryland (College Park).

# <a name="InstallGuide"></a>Getting Started

This section provides a very brief guide to get started with PySedSim. Much more detailed instructions appear in the PySedSim user manual, especially regarding the preparation of input data to the model. The user manual is available in the /docs directory of this repository. Users may wish to use the PySedSim example cases provided in the /example directory to test the steps below, as well as to see how to create the input model files for your own model application.

1.	Download and install Python 2. (Note: v. 2.7 or later is suggested for Python 2; PySedSim has not been tested with Python 3). 

You can get a free distribution of Python and most of its most popular packages (libraries) for scientific computing through Anaconda: https://www.continuum.io/downloads.

2. Create a directory on your workstation into which to clone PySedSim. We will refer to this as pysedsim_main.

3. We recommend conducting the clone process (discussed in the following steps) from the command line. If you are unfamiliar working from the command line, we suggest you downloading Git for free here: https://git-scm.com/downloads. This Git download will come with “git bash”, from which you can attempt the “git clone” operation described below. 

4. This repository uses the Git Large File Storage (LFS) extension (see https://git-lfs.github.com/ for details). Please install GitLFS and run the following command in your pysedsims_main directory before cloning if you do not already have Git LFS initialized: 

>>git lfs install.

5. Clone PySedSim into your pysedsim_main directory by entering the following command in a prompt (e.g., Git bash):

>> git clone https://github.com/FeralFlows/PySedSim.git

The above command should initiate a download of PySedSim files, including large files stored with Git LFS. If none of the large files download in the example/ directory, windows users have have better luck with:

>> git lfs clone https://github.com/FeralFlows/PySedSim.git

6. From the directory you cloned PySedSim into, run:

>>python setup.py install . 

This will install PySedSim as a Python package on your machine and install all of the needed dependencies. If you receive any errors that certain packages are not installed, we recommend installing these packages through the Anaconda command prompt by typing the command below (inserting the correct package where “package_name” appears”).

>>pip install package_name

7. Run the provided example simulation(s).

We recommend first running the example formulations provided out-of-the-box in the cloned repository. To do this, change directories to a formulation’s directory (e.g., formulation 1):

>>cd pysedsim_sims_main/example/formulation_1

Then, run pysedsim from the command prompt with the following command:

>>python formuation_1.py

Check the output directories, and log file, to confirm that the simulation has run correctly.

8. Setup and run your own simulation.

We refer you to the user manual in this repository's /docs directory for images and examples of how to prepare and correctly format input files for your own PySedSim application.

# <a name="ModelFileStructure"></a>Overview of Model File Structure

The figure below summarizes the three steps involved in running and reviewing a PySedSim simulation: preparing input files (upper), running a simulation (middle), and viewing outputs (lower). As shown in Figure 4 1 (upper), to run a simulation (or simulation-optimization) experiment the user is required specify at least two input files: (1) a top-level “configuration.csv” input file, referred to henceforth as the configuration file; and (2) a main “Input_Data_File.xlsx” input file, referred to henceforth as the main input data file. In the event the user is running a Monte Carlo simulation, Figure 4 1 (upper) shows that the user may also specify two file types: 1) “.csv” files containing time series of stochastic realizations of daily hydrologic and sediment inflows at each incremental inflow junction (see description in Appendix A of user manual), with a unique file for every junction, and water and sediment appearing in different files; and 2) a “Monte_Carlo_Input_Specifications.xlsx”, referred to henceforth as the Monte Carlo input file. The configuration file is used to specify the file paths of all input and output files, as well as to run multiple experiments in batch mode. The main input data file serves as the model’s user interface, giving the user control over model operation as well as data input, editing, and output, and in doing so facilitates PySedSim’s use as a decision support system. This file organizes related inputs into worksheets (i.e., tabs).

![pysedsim_schematic](/images/pysedsim_schematic.png)

# <a name="Pubs"></a>Publications using PySedSim

<strong> Peer-reviewed Publications: </strong>

Wild, T.B., Reed, P.M., Loucks, D.P., Mallen-Cooper, M., Jensen, E.D. (2019). Balancing Hydropower and Ecological Impacts in the Mekong: Tradeoffs for Sambor Mega Dam. J. Water Resour. Plann. Manage. DOI: 10.1061/(ASCE)WR.1943-5452.0001036.

<strong> Presentations: </strong>

Loucks, D.P., Wild, T.B., Reed, P.M (2018). Ecology-Energy Tradeoffs in the Lower Mekong. Oral presentation at International Symposium on Safety of Water Resources Engineering for Regulation of water Cycle in River Basin and Disaster Prevention and Mitigation, IWHR, Beijing, China, October 2018.

Loucks, D.P., Wild, T.B., Reed, P.M (2018). Communicating complex hydropower-ecological tradeoffs and their uncertainties to decision makers: A hydroinformatics challenge in the Lower Mekong. Oral presentation at European Geosciences Union (EGU) General Assembly 2018, Vienna, Austria. 

Loucks, D.P., Wild, T.B., and Reed, P.M. (2017). Tradeoffs among Hydropower, Sediment Flow, and Fish Survival in the Lower Mekong: A Study of Sambor Dam. Hydropower Sustainability Forum: Mekong +, Multiconsult, Oslo, Norway, September 2017.

Loucks, D.P., Wild, T.B., and Reed, P.M. (2017). Modeling and Managing Water in a Changing and Surprising World. Oral presentation at ASCE World Environmental and Water Resources Congress 2017, Sacramento, CA.

Wild, T.B., Reed, P.M., and Loucks, D.P. (2016). Mega Dams on the Mekong: A Framework for Evaluating Alternative Dam Siting, Design and Operations Options to Improve Sediment and Migratory Fish Passage. Oral presentation at ASCE World Environmental and Water Resources Congress 2016, West Palm Beach, FL.

Wild, T.B., Reed, P.M., and Loucks, D.P. (2015). Food-Energy Tradeoffs in the Mekong: Options for Improved Sediment and Fish Passage at Hydropower Dams. Oral presentation at ASCE World Environmental and Water Resources Congress 2015, Austin, TX.

Loucks, D.P. and Wild, T.B. (2015). Managing Tradeoffs between Hydropower and the Environment in the Mekong River Basin. Oral presentation at European Geosciences Union (EGU) General Assembly 2015, Vienna, Austria.
