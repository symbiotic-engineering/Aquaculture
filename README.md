# Aquaculture
This repository contains an "optimal design of offshore aquaculture with co-located wave energy", simulated through Python.

## File Structure
This repository consists of the following files:
- modules: scripts and functions to define "Wave-powered Aquaculture" optimization problem.
- objects: scripts and functions to define aquaculture and wave energy converter objects and their parameters.
- optimization: scripts and functions to perform a single objective optimization.
- wpaf_opt: top-level script to perform the "WPAF" optimization code.
- run_sim_wpaf: top-level script to set the required input files for wpaf_opt.
- gis/gis_handler: object to handle connection between GIS data and Python functions
- gis/gis_example: examples of how to use the GIS handler to query and save data
- gis/data: GIS files used in analysis including conflicts (.geojson) and conditions (.tif)
- requirements.txt: list of the required libraries.

## How to use
In order to run the optimization code and find optimal design for wave-powered aquaculture farm (WPAF), you can open and run the run_sim_wpaf (.ipynb) file.

## Dependencies
Python 3.8.8 is used to develop this package. 
For other required libraries, please refer to requirements.txt file.

## Context
The project is part of research in the Symbiotic Engineering Analysis (SEA) Lab and will be submitted to the Renewbale Energy.

## Authors
- Arezoo Hasankhani, ah844@cornell.edu (Point of contact)
- Gabriel Ewig, gre27@cornell.edu
- Eugene Thome Won, etw36@cornell.edu
- Maha Haji, maha@cornell.edu

## Funding Acknowledgement
This material is based upon work supported by the Sea Grant Regional Research Project No.: R/ATD-18-NESG.

## License
This project is released open-source under the MIT License. Processed GIS data included in the repository is from publicly available sources, but may have different license terms.
