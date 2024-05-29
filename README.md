# TumInvasion-BM: Simulation of the rupture of the basement membrane by the effect of tumor cells

## Description
This repository contains the codes for the simulations of the rupture of the basement membrane (BM).
BM rupture happens because of the effect of digestive enzymes activated by the tumor cells.
These active enzymes are produced in inactive form by the senescent fibroblasts (SFs).

The structure of the repository is the following:
- parameters.m: contains the parameters values.
- reduced_system/: contains codes for the simulations of the reduced system.
  + fun_parameters.m: contains the definitions of the important parameters of the model.
  + fun_Sim_ReducedModel.m: contains the functions to simulate the PDE model and the ODE system.
  + Main_Plot.m: is the main script. 

## How to use
### Simulating the model
For the reduced and full systems, the code is run from the main scripts "script_reduced_system.m" and "script_full_system.m" respectively.
In the beginning of these scripts, the user must choose:
- the dimension of the conjunctive tissue dim = 1 or 2. The BM will have the dimension dim-1. 
- the test: if test = 1 it corresponds to the healthy test case (no tumor cells nor SF), 
test = 2 corresponds to the tumor test case (tumor cells in the BM and no SF), 
test = 3 corresponds to the tumor + SF case (positive tumor cells density and SF in the conjunctive)

The scripts are composed of the following sections: 
- discretization (space and time)
- parameter values (in this section there is a call to the parameters.m function). 
After this call, there is the part to prepare the source terms and to include the effect of the SFs in the model. 
- initial conditions. Here is a call to the initial_conditions.m file. It prepares the vectors corresponding the initial concentrations of enzymes and densities.
- Stability test for the reduced system. It tests if the stability of the ODE system (see article's supplementary information)
- Call the main function to solve the system. 

### Plotting
If you want to obtain plots of the quantities, run the plot_figs.m script after the completion of the main script. 
At the beginning of the plot_figs.m script, you can select if you want to save the figures (to the directory specified in "plotpath" variable)

## Authors
Alexandre Poulain and Chiara Villa. 

## How to cite
Please cite our article if you are using this code: Luís Almeida, Alexandre Poulain, Albin Dr Pourtier, Chiara Villa. Mathematical modelling of the contribution of senescent fibroblasts to basement membrane digestion during carcinoma invasion. 2024. ⟨hal-04574340⟩
