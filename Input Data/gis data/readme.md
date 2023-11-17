# Data Files Folder

This folder contains all necessary data files for the Python handler and model including rasters of oceanic conditions, vectors of potential conflicts, and a vector file with the included state and federal waters.

Due to space limitations on GitHub, these files are hosted externally on Zenodo and must be downloaded separately. The current files can be downloaded from [Zenodo](https://zenodo.org/records/10140826) but may have minor deviations from those initially used in published results due to additional pre-processing and scaling. Once downloaded, the files should be placed in this directory with the following paths:

    ### Conditions
    /Input Data/gis data/Surface Currents m-s (NODP 2016).tif
    /Input Data/gis data/Surface Oxygen mg-l (NCEI 2019).tif
    /Input Data/gis data/Surface Salinity PSU (NCEI 2019).tif
    /Input Data/gis data/Surface Temperature C (NODP 2016).tif
    /Input Data/gis data/Wave Energy Period s (NREL 2011).tif
    /Input Data/gis data/Significant Wave Height m (NREL 2011).tif
    /Input Data/gis data/Bathymetry Downsampled m (NGDC 1990).tif
    /Input Data/gis data/Distance to Port m (OCM 2019).tif
    /Input Data/gis data/Distance to Shore m (OCM 2018).tif

    ### Conflicts
    /Input Data/gis data/Very High Fishing Vessel Traffic (NODP 2022).geojson
    /Input Data/gis data/High Fishing Vessel Traffic (NODP 2022).geojson
    /Input Data/gis data/Marine Protected Areas (NMPAC 2020).geojson
    /Input Data/gis data/Danger Zones and Restricted Areas (OCM 2022).geojson
    /Input Data/gis data/Submarine Transit Lanes (NODP 2016).geojson
    /Input Data/gis data/Cape Cod TORPEX (NODP 2016).geojson
    /Input Data/gis data/Block Island Renewable Energy Zone (NODP 2010).geojson
    /Input Data/gis data/MA Wind Energy Areas (NODP 2015).geojson
    /Input Data/gis data/Potential Wind Lease Areas (BOEM 2023).geojson
    /Input Data/gis data/Wind Planning Areas (BOEM 2023).geojson
    /Input Data/gis data/Shipping Lanes (OCS 2015).geojson

    ### Scope
    /Input Data/gis data/Northeast State and Federal Waters (OCM 2018).geojson