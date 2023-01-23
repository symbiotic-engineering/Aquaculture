# Aquaculture
This repository contains an "Optimization of offshore aquaculture with co-located wave energy", simulated through Python.

## File Structure
This repository consists of the following files:
- modules: scripts and functions to define "Wave-powered Aquaculture" optimization problem.
- objects: scripts and functions to define aquaculture objects and their parameters.
- optimization: scripts and functions to perform a single objective optimization using SciPy open-source library.
- run_sim_wec&pen: top-level script to perform the "Wave-powered Aquaculture" optimization code using default values for environmental parameters.
- run_sim_wec&pen_random_init: top-level script to perform the "Wave-powered Aquaculture" optimization code using random values for starting points of design variables to check the convergence of the optimization.
- param_sweep: top-level script to perform the "sensitivity analysis".

## How to use
Sensitivity Analysis: In order to run the "sensitivity analysis", you can open and run the param_sweep (.ipynb) file.

In order to run the "Wave-powered Aquaculture" optimization code and find optimal wave energy convereter (WEC) and optimal pen considering constant environmental parameters, you can open and run the run_sim_wec&pen (.ipynb) file.

In order to run the "Wave-powered Aquaculture" optimization code and find optimal WEC and optimal pen with random values for starting points of design variables to check the convergence of the optimization, you can open and run the run_sim_wec&pen_random_init (.ipynb) file.

## Dependencies
Python 3.8.8 is used to develop this package. 
For other required libraries, please refer to requirements.txt file.

## Context
The project is part of research in the Symbiotic Engineering Analysis (SEA) Lab and has been submitted to the International Society of Offshore and Polar Engineers (ISOPE 2023).

## Authors
- Arezoo Hasankhani, ah844@cornell.edu (Point of contact)
- Rebecca McCabe, rgm222@cornell.edu 
- Gabriel Ewig, re27@cornell.edu
- Eugene Thome Won, etw36@cornell.edu
- Maha Haji, maha@cornell.edu

## Funding Acknowledgement
This material is based upon work supported by the Sea Grant Regional Research Project No.: R/ATD-18-NESG.

## License
This project is released open-source under the MIT License. 
