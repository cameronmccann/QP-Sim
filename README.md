# Evaluating Propensity Score Models for Estimating Causal Mediation Effects under Unmeasured Cluster-level Confounding

## Overview

We conducted two simulation studies and an empirical application to assess the performance of three propensity score (PS) models in estimating mediation effects with clustered data when an unmeasured cluster-level confounder is present. This repository contains the code and data used at each stage of the study (i.e., Simulation Study 1, Simulation Study 2, and Empirical Application).

In Simulation Study 1 and 2, there are individual-level treatment (binary), mediator (continuous), and outcome (continuous) variables. In Simulation Study 1, we assessed the PS models’ ability to accurately estimate the indirect effect when both an unmeasured cluster-level confounder and three individual-level covariates influence the treatment-mediator and mediator-outcome relationships. 

In Simulation Study 2, we modified the simulation design so that the unmeasured cluster-level confounder influences only the treatment-outcome relationship. The three individual-level covariates still influence both the treatment-mediator and mediator-outcome relationships; however, an additional set of three covariates was introduced to influence only the treatment-outcome relationship. Simulation Study 2 also includes two versions of the direct and indirect effect estimates: the Total Natural Direct Effect (TNDE), the Pure Natural Direct Effect (PNDE), the Total Natural Indirect Effect (TNIE), and the Pure Natural Indirect Effect (PNIE).

The Empirical Application uses data from the National Longitudinal Study of Adolescent to Adult Health (Add Health; Harris & Udry, 2008) to illustrate the application of the methods evaluated. The cleaned version of the dataset used is included in the `Application/Data/Cleaned` folder.


## Folder Structure

The repository is organized into four main folders, which include the following information:
```
├── README.md    
├── Code/                          # Scripts to conduct simulations (i.e., generating & analyzing data) & obtain results (i.e., RMSE & relative bias)
│   ├── S1_Conduct-Simulation.R   
│   ├── S1_Obtain-Results.R       
│   ├── S2_Conduct-Simulation.R   
│   ├── S2_Obtain-Results.R       
├── Functions/                     # Helper functions for simulations (generating data & analyzing data)  
│   ├── AnalysisFunc_Sim1.R
│   ├── genOneData_Sim1.R
│   ├── ...
├── Output/                        
│   ├── S1_Results/                # Performance measures for each simulation condition 
│   │   ├── Data/
│   │   └── Figures/ 
│   ├── S1_Simulation-Output/      # Simulation output (e.g., direct & indirect estimates) for each run under each simulation condition
│   ├── ...
└── Application/                   # Main folder for the empirical application 
    ├── Code/                      # Scripts for data cleaning and analysis
    ├── Data/                      # Raw and cleaned datasets
    ├── Functions/                 # Functions to perform bootstrap confidence intervals 
    └── Output/                    # Empirical application results, visualizations, estimates, and temporary boostrap data files
        ├── Visuals/               
        ├── Estimates/
        └── Bootstrap_Temp/        
```

## Requirements and Dependencies

This repository requires R (version 4.x or higher) and specific R packages for conducting the simulation studies and the empirical application. The scripts use the pacman package to automatically install and load missing packages, but you can manually install them if preferred. 

### Simulation Study
#### Conducting Simulations (1 and 2)

The following packages are required to run the simulations:

- `doParallel`
- `foreach`
- `parallel`

#### Obtain Simulation Results (1 and 2)

The following packages are required to analyze the simulation results:

- `tidyverse`
- `ggplot2`
- `flextable`
- `stringr`
- `ggdag`
- `dagitty`
- `huxtable`

### Empirical Application

The following packages are required for the empirical application:

- `tidyverse`
- `ggplot2`
- `extrafont`
- `stringr`
- `mice`
- `cobalt`
- `WeightIt`
- `boot`
- `utils`
- `lme4`
- `WeMix`
- `parallel`


## How to Use This Repository

### Running Simulations

To conduct the simulations and obtain results, follow these steps: 

1. Clone or download the repository.
2. Navigate to the `Code` folder.
3. Run the simulation scripts in the following order:
    - `S1_Conduct-Simulation.R`: Runs the first simulation and stores output (e.g., direct & indirect estimates) for each run under each simulation condition in the `Output/S1_Simulation-Output` folder.
    - `S1_Obtain-Results.R`: Processes output from the first simulation by computing performance metrics (e.g., relative bias) and generating visuals (stored in `Output/S1_Results`).
    - Repeat the process for the second simulation (`S2`).
        
        
### Running the Empirical Application

1. Navigate to the `Application/Code` folder.
2. Run `Empirical-Application_Data-Cleaning.R` to prepare the data for analysis.
3. Execute `Empirical-Application_Analysis.R` to analyze the empirical dataset and obtain the effect estimates.
4. The results, including visualizations and estimates, will be saved in the `Application/Output` folder.

        
