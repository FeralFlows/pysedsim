# PySedSim
An Open Source Reservoir and Sediment Simulation Screening Framework for Identifying and Evaluating Dam Siting, Design, and Operation Alternatives

- [Contact Us](#Contact)
- [Introduction](#Introduction)
- [Getting Started](#InstallGuide)
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

1. Clone SedSim from Github.

   The installation process is described in detail in the SedSim user manual in the docs folder. The basic steps are as follows:

   Users can either directly download SedSim files from the model’s github repository (https://github.com/FeralFlows/SedSim), or can “git clone” SedSim from the command prompt with the following command:

   “git clone https://github.com/FeralFlows/SedSim.git”

   If you are unfamiliar working from the command line, we suggest you try downloading Git for free here: https://git-scm.com/downloads. This download will come with “git bash”, from which you can attempt the git clone operation mentioned above. Cloning files rather than downloading them directly will enable you to pull recent SedSim repository updates directly to the local repository on your computer.

2. Open up the main (macro) workbook (SedSim.xlsm).

   Open the main SedSim model file. This file will be referred to as "SedSim_Model.xlsm" throughout this documentation for convenience, but the file can be given any name by simply right-clicking on the file icon when the file is closed and selecting "rename".
Upon opening the workbooks, you may be asked if you wish to enable macros. Click the “Enable Macros” option to allow the sediment model to execute properly.

   The image below shows the interface associated with the main workbook of the SedSim model.  It is designed to be generic so it can be used without code modification to run any input file. 

![sedsim_worksheet](/images/Capture.PNG)

3. Enable macros in security settings.

   To be certain that the SedSim model will always run on your computer, in the “SedSim_Model.xlsm” workbook, in Excel 2007 (or Excel 2010), go to File (or MS Office Button) -> Options -> Trust Center -> Trust Center Settings -> Macro Settings -> Enable All Macros. When you are finished running the model in Excel, these settings should be returned to their original status to avoid potential security threats to your computer. Alternatively, as described in step 1 above, your version of Excel may provide a warning message when you first open the SedSim model that asks if you wish to enable the currently open SedSim model Excel file to be run on your computer, among other Macro options. You can enable the model to be run on your computer this way as well. 

   The next two figures below visually depict the steps described above. From the options menu within Excel’s “File” tab, select “Trust Center”, from which you can enable macros.
   
   ![trust_center](/images/trust_center.png)
   ![trust_center](/images/trust_center_2.png)

4. Load in the input data and specify assumptions.

   Load your simulation assumptions and data into the main input file (“SedSim_Input.xlsx”) files. You can obtain an empty version of this file on the SedSim Github repository. Alternatively, you can copy the input file from the “example” directory on the paper’s Github repository (SedSim_Input_Example.xlsx) and use it as an example, replacing the example data with your data. All colored worksheets will require some input, whereas un-colored worksheets will not require user input and are instead populated during the execution of the macro. Model-related assumptions (e.g., sediment density) can be modified in the "Simulation Specifications" worksheet of the input data workbook. Please review the “Input File” section of this user manual for specific details regarding how to populate each worksheet within the main input file.

   The “SedSim_Model.xlsm” workbook is the only file that is required to be open for the simulation to run properly. The Input file and output file do not need to be opened beforehand. The input file will be opened and closed automatically when needed by the main macro, and the output file is automatically created, saved, and closed by the main macro.

5. Run the model. 

   This can be performed by clicking the "Run Model" button in the "SedSim_Model.xlsm" workbook.  During the execution of the model, the model will automatically close the input data file. The model was designed to automatically close the input data file once all data have been imported into internal arrays because keeping the input file open can exhaust the maximum memory usage limits of Excel for a large reach/reservoir network or long simulation duration. The model may produce two different types of error messages during execution: (1) a detailed error message generated by SedSim that the user must acknowledge, by clicking "OK" on the automatically generated error message box, before the simulation can proceed; and (2) an excel VBA error message, which is not likely to contain detailed instructions, and which is likely the result of improper input data specification or input/output file naming.

   Note: If the model will not run and displays an error regarding the Microsoft Excel “Solver” package, you may need to install “References” within VBA, as described in Step 5.

6. If you experience errors during a model run, install necessary “References” within Excel Visual Basic.

   Do this by opening up the "SedSim_Model.xlsm" and accessing the VBA code by selecting Alt+F11 on the keyboard. Within the “Project” menu on the left-hand side of the screen, select (VBAProject (SedSim.xlsm), click on the main model file to reveal its sub-menu, then double click the “SedSim_Model” module within the sub-menu. 

   In the main menu at the top of the screen, select Tools-->References. Find and check two boxes: (1) “Solver”, which enables sediment calibration; and (2) “Microsoft Scripting Runtime”, which enables runtime messages to be printed to a text file during model execution. Click OK to install the solver references. This is shown in the two images below.
 
   After installing these references, save the file within VBA, and close out of VBA and Excel altogether. Finally, re-open and re-run the model to see if this change permits the model to run.

   ![trust_center](/images/references.png)
   ![trust_center](/images/references_2.png)

7. Evaluate results.   

   The results of the simulation run are contained in the “SedSim_Output.xlsx” file.  Are these results reasonable given the input data?  One approach to gain confidence in the results is to create input data for relatively simple systems that should lead to obvious results, and then see if indeed they did.

# <a name="Pubs"></a>Publications using PySedSim

<strong> Peer-reviewed Publications: </strong>

Wild, T.B., Reed, P.M., Loucks, D.P., Mallen-Cooper, M., Jensen, E.D. (2018). Balancing Hydropower and Ecological Impacts in the Mekong: Tradeoffs for Sambor Mega Dam. J. Water Resour. Plann. Manage. DOI: 10.1061/(ASCE)WR.1943-5452.0001036.

<strong> Presentations: </strong>

Loucks, D.P., Wild, T.B., Reed, P.M (2018). Ecology-Energy Tradeoffs in the Lower Mekong. Oral presentation at International Symposium on Safety of Water Resources Engineering for Regulation of water Cycle in River Basin and Disaster Prevention and Mitigation, IWHR, Beijing, China, October 2018.

Loucks, D.P., Wild, T.B., Reed, P.M (2018). Communicating complex hydropower-ecological tradeoffs and their uncertainties to decision makers: A hydroinformatics challenge in the Lower Mekong. Oral presentation at European Geosciences Union (EGU) General Assembly 2018, Vienna, Austria. 

Loucks, D.P., Wild, T.B., and Reed, P.M. (2017). Tradeoffs among Hydropower, Sediment Flow, and Fish Survival in the Lower Mekong: A Study of Sambor Dam. Hydropower Sustainability Forum: Mekong +, Multiconsult, Oslo, Norway, September 2017.

Loucks, D.P., Wild, T.B., and Reed, P.M. (2017). Modeling and Managing Water in a Changing and Surprising World. Oral presentation at ASCE World Environmental and Water Resources Congress 2017, Sacramento, CA.

Wild, T.B., Reed, P.M., and Loucks, D.P. (2016). Mega Dams on the Mekong: A Framework for Evaluating Alternative Dam Siting, Design and Operations Options to Improve Sediment and Migratory Fish Passage. Oral presentation at ASCE World Environmental and Water Resources Congress 2016, West Palm Beach, FL.

Wild, T.B., Reed, P.M., and Loucks, D.P. (2015). Food-Energy Tradeoffs in the Mekong: Options for Improved Sediment and Fish Passage at Hydropower Dams. Oral presentation at ASCE World Environmental and Water Resources Congress 2015, Austin, TX.

Loucks, D.P. and Wild, T.B. (2015). Managing Tradeoffs between Hydropower and the Environment in the Mekong River Basin. Oral presentation at European Geosciences Union (EGU) General Assembly 2015, Vienna, Austria.
