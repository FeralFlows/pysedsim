# PySedSim
An Open Source Reservoir and Sediment Simulation Screening Framework for Identifying and Evaluating Dam Siting, Design, and Operation Alternatives

![sedsim_worksheet](/images/pysedsim_logo.png)

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
