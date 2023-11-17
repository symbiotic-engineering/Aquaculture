# Aquaculture
This repository contains an "optimal design of offshore aquaculture with co-located wave energy and energy storage", simulated through Python.

## File Structure
This repository consists of the following files:
- modules: scripts and functions to define "Wave-powered Aquaculture" optimization problem.
- objects: scripts and functions to define aquaculture, wave energy converter, and energy storage objects and their parameters.
- optimization: scripts and functions to perform a single objective optimization, and multi-objective optimization.
- wpaf_opt: top-level script to perform the "WPAF" optimization code and set the required input files.
- run_sim_soo_wpaf: top-level script to run the single objective optimization.
- run_sim_moo_wpaf: top-level script to run the multi-objective optimization.
- run_sim_soo_wpaf_random_init: top-level script to run the single objective optimization with random initial points for the design variables to check the optimization convergence.
- sensitivity_analysis: script to perform sensitivity analysis of the objective to the design variables.
- utilities: script to print and draw figures.
- requirements.txt: list of the required libraries.

## How to use
- In order to run the single objective optimization code and find the optimal design for a wave-powered aquaculture farm (WPAF), you can open and run the run_sim_soo_wpaf (.ipynb) file.
- In order to run the single objective optimization code with random initial points for the design variables and check the convergence of the optimization, you can open and run the run_sim_soo_wpaf_random_init (.ipynb) file.
- In order to run the multi-objective optimization code and find the optimal design for a wave-powered aquaculture farm (WPAF), you can open and run the run_sim_moo_wpaf (.ipynb) file.
- In order to run the sensitivity analysis and find the sensitivity of the objective function to the design variables, you can open and run the sensitivity_analysis (.ipynb) file.


## Dependencies
Python 3.8.8 is used to develop this package. 
For other required libraries, please refer to requirements.txt file.

## Context
The project is part of research in the Symbiotic Engineering Analysis (SEA) Lab and will be submitted to the Renewbale Energy.

## Citation
Hasankhani, A., McCabe, R., Ewig, G., Won, E. T., & Haji, M. N. (2023, June). Conceptual design and optimization of a wave-powered offshore aquaculture farm. In The 33rd International Ocean and Polar Engineering Conference. OnePetro.

## Authors
- Arezoo Hasankhani, ah844@cornell.edu (Point of contact)
- Eugene Thome Won, etw36@cornell.edu
- Maha Haji, maha@cornell.edu

## Funding Acknowledgement
This material is based upon work supported by the Sea Grant Regional Research Project No.: R/ATD-18-NESG.

## License
This project is released open-source under the MIT License. Processed GIS data included in the repository is from publicly available sources, but may have different license terms.
