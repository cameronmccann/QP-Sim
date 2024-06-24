# Causal Mediation with Unmeasured Cluster-level Confounders: Evaluating Propensity Score Models


## Overview

This repository contains the code and documentation for a simulation study on [brief description of the study topic]. The study aims to [main objective of the study].

## Table of Contents

## Folder Structure

The main 3 folders provide the following information:
- `Code` contains scripts to run simulations and report results (i.e., compute and display performance measures).
- `Functions` contains functions that are implemented when running simulations (e.g., generate data).
- `Output` contains subfolders holding simulation output information (e.g., estimates per run) & results (e.g., performance measures). 

The repository is organized as follows:
```
├── QP-Sim.Rproj
├── README.md    
└── Code
    ├── S1_Conduct-Simulation.R   # Simulation 1 code (i.e., generating & analyzing data)
    ├── S1_Obtain-Results.R       # Script for reporting RMSE & relative bias
    ├── S2B_Conduct-Simulation.R  # Simulation 2B code (i.e., generating & analyzing data)
    ├── S2B_Obtain-Results.R      # Script for reporting RMSE & relative bias
    ├── S2_Conduct-Simulation.R   # Simulation 2 code (i.e., generating & analyzing data)
    ├── S2_Obtain-Results.R       # Script for reporting RMSE & relative bias
└── Functions
    ├── AnalysisFunc_Sim1.R
    ├── AnalysisFunc_Sim2.R
    ├── AnalysisFunc_Sim2B.R
    ├── genOneData_Sim1.R
    ├── genOneData_Sim2.R
    ├── genOneData_Sim2B.R
└── Output
    └── S1_Results/Data
    └── S1_Simulation-Output
    └── S2B_Results/Data
    └── S2B_Simulation-Output
    └── S2_Results/Data
    └── S2_Simulation-Output
```



