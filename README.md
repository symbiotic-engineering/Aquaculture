[![DOI](https://zenodo.org/badge/467335021.svg)](https://zenodo.org/doi/10.5281/zenodo.7633737)
![GitHub](https://img.shields.io/github/license/symbiotic-engineering/aquaculture)

# Aquaculture
This repository contains a Python to model an offshore aquaculture farm with co-located wave energy converters and determine potential sites using marine spatial planning.

## File Structure
This repository consists of the following files:
- modules: scripts and functions to define "Wave-powered Aquaculture" optimization problem.
- objects: scripts and functions to define aquaculture and wave energy converter objects and their parameters.
- optimization: scripts and functions to perform a single objective optimization using brute-force.
- env_bruteforce: script to perform brute-force for defined net pen parmeters.
- run_sim_env: top-level script to perform the "marine spatial planning" optimization code using default values for defined environment with rasters and vectors of data.
- gis_handler: object to handle connection between GIS data and Python functions
- scripts/: relevant scripts including examples of how to use the GIS handler to query and save data
- data/: GIS files used in analysis including conflicts (.geojson) and conditions (.tif) *(stored externally on Zenodo)*
- requirements.txt: list of the required libraries.

## How to use
Relevant data files must be first downloaded from the Zenodo at [https://zenodo.org/records/10140826](https://zenodo.org/records/10140826) and placed in the data folder of this repository. In order to run the marine spatial planning optimization code and find optimal location for wave-powered aquaculture farm (WPAF) for constant net pen within the Northeast U.S, open and run the run_sim_env (.ipynb) file.

## Dependencies
Python 3.9.18 is used to develop this package.
For other required libraries, please refer to requirements.txt file.

## Context
The project is part of research in the Symbiotic Engineering Analysis (SEA) Lab and is published in the Oceans 2023.

## Citation
Hasankhani, A., Ewig, G., McCabe, R., Won, E. T., & Haji, M. (2023, June). Marine Spatial Planning of a Wave-Powered Offshore Aquaculture Farm in the Northeast US. In OCEANS 2023-Limerick (pp. 1-10). IEEE.

## Authors
- Arezoo Hasankhani, ah844@cornell.edu (Point of contact)
- Gabriel Ewig, gre27@cornell.edu
- Rebecca McCabe, rgm222@cornell.edu 
- Eugene Thome Won, etw36@cornell.edu
- Maha Haji, maha@cornell.edu

## Funding Acknowledgement
This material is based upon work supported by the Sea Grant Regional Research Project No.: R/ATD-18-NESG.

## License
This project is released open-source under the MIT License. Processed GIS data included in the repository is from publicly available sources, but may have different license terms.
