# Causal Mediation with Unmeasured Cluster-level Confounders: Evaluating Propensity Score Models


## Overview

This repository contains the code and documentation for the present simulation study. The study aims to assess the performance of three different propensity score (PS) models for estimating mediation effects in mediation analysis with multilevel data when treatment assignment is nonrandomized and an unmeasured cluster-level confounder exists. 

## Folder Structure

The main 3 folders provide the following information:
- `Code` contains scripts to run simulations and report results (i.e., compute and display performance measures).
- `Functions` contains functions that are implemented when running simulations (e.g., generate data).
- `Output` contains subfolders holding simulation output information (e.g., estimates per run) & results (e.g., performance measures). 

The repository is organized as follows:
```
├── QP-Sim.Rproj
├── README.md    
└── Code                          # Code to conduct simulations (i.e., generating & analyzing data) & report results (i.e., RMSE & relative bias)
    ├── S1_Conduct-Simulation.R   
    ├── S1_Obtain-Results.R       
    ├── S2B_Conduct-Simulation.R  
    ├── S2B_Obtain-Results.R      
    ├── S2_Conduct-Simulation.R   
    ├── S2_Obtain-Results.R       
└── Functions                     # Functions used in simulations 
    ├── AnalysisFunc_Sim1.R
    ├── AnalysisFunc_Sim2.R
    ├── AnalysisFunc_Sim2B.R
    ├── genOneData_Sim1.R
    ├── genOneData_Sim2.R
    ├── genOneData_Sim2B.R
└── Output                        # Output from simulations 
    └── S1_Results/Data
    └── S1_Simulation-Output
    └── S2B_Results/Data
    └── S2B_Simulation-Output
    └── S2_Results/Data
    └── S2_Simulation-Output
```

The repository is organized as follows:
```
├── QP-Sim.Rproj
├── README.md    
├── Code                          # Code to conduct simulations (i.e., generating & analyzing data) & report results (i.e., RMSE & relative bias)
│   ├── S1_Conduct-Simulation.R   
│   ├── S1_Obtain-Results.R       
│   ├── S2B_Conduct-Simulation.R  
│   ├── S2B_Obtain-Results.R      
│   ├── S2_Conduct-Simulation.R   
│   ├── S2_Obtain-Results.R       
├── Functions                     # Functions used in simulations 
│   ├── AnalysisFunc_Sim1.R
│   ├── AnalysisFunc_Sim2.R
│   ├── AnalysisFunc_Sim2B.R
│   ├── genOneData_Sim1.R
│   ├── genOneData_Sim2.R
│   ├── genOneData_Sim2B.R
├── Output                        # Output from simulations 
│   ├── S1_Results
│   │   └── Data
│   ├── S1_Simulation-Output
│   ├── S2B_Results
│   │   └── Data
│   ├── S2B_Simulation-Output
│   ├── S2_Results
│   │   └── Data
│   ├── S2_Simulation-Output
│   └── Estimates
│       └── Estimates.rds
└── Application                   # Main folder for the empirical application 
    ├── Code                      # Scripts for data cleaning and analysis
    │   ├── Empirical-Application_Analysis.R
    │   └── Empirical-Application_Data-Cleaning.R
    ├── Data                      # Folder containing datasets
    │   ├── Raw                   # Original datasets (no modifications)
    │   │   └── ICPSR_21600       # Subfolder for specific original data sources
    │   └── Cleaned               # Final clean dataset used for analysis
    │       └── Empirical-Application-Data.rds
    ├── Functions   
    │   ├── bootstrap_ci_paral_2.R
    │   ├── bootstrap_ci_paral.R
    │   ├── bootstrap_ci_re_mean_paral.R
    │   ├── bootstrap_ci_re_paral_2.R
    │   └── bootstrap_ci_re_paral.R
    └── Output                    # Results from the empirical analysis
        ├── Visuals               # Plots, graphs, and other visuals
        ├── Estimates
        │   └── Effect-Estimates.rds
        └── Bootstrap_Temp        # Temporary bootstrap data files
```


