# mtDesertIsland_Pinusrigida
Data and scripts corresponding to Pinus rigida study from Mount Desert Island, Maine, USA

## general information
This repository houses data and scripts from a study in Mount Desert Island, Maine, USA examining Pinus rigida traits in high vs. low elevation and fire vs no fire history. 

## authors
Jeff Licht and Nick Smith (nick.gregory.smith@gmail.com)

## DOI badge
[![DOI]()

## folder descriptions
### /data
This folder contains all of the data files.

### /scripts
This folder contains all of the R scripts.

## file descriptions
### /data/mdi_all_clean.csv
This file contains topography, allometry, soil, and foliar data. The specific column information is:
- *ID*: the sample ID indicating site_rep_elevation
- *Name*: site name
- *longitude*: longitude (decimal degrees) at individual tree
- *latitude*: latitude (decimal degrees) at individual tree
- *Elevation*: elevation (m) at individual tree
- *Slope*: slope at individual tree measured using ArcGIS 
- *Aspect*: aspect at individual tree measured using ArcGIS
- *Height*: tree height measured with nested, 2 m calibrated, aluminum rods (m)
- *Canopy*: canopy spread measured with nested, 2 m calibrated, aluminum rods (m)
- *Diam*: diameter at breast height measured with a ProSkit electronic digital caliper (cm)
- *d13C*: foliar carbon-13 isotope measured with a Thermo Delta V+ IR-MS continuous flow isotope ratio mass spectrometer (‰)
- *d15N*: foliar nitrogen-15 isotope measured with a Thermo Delta V+ IR-MS continuous flow isotope ratio mass spectrometer (‰)
- *C_foliar*: foliar carbon measured with a Leco CN-2000 Carbon-Nitrogen Analyzer (g<sup>1</sup> g<sup>-1</sup>)
- *N_foliar*: foliar nitrogen measured with a Leco CN-2000 Carbon-Nitrogen Analyzer (g<sup>1</sup> g<sup>-1</sup>)
- *Ca_foliar*: foliar calcium (g<sup>1</sup> g<sup>-1</sup>)
- *P_foliar*: foliar phosphorus (g<sup>1</sup> g<sup>-1</sup>)
- *K_foliar*: foliar potassium (g<sup>1</sup> g<sup>-1</sup>)
- *Mg_foliar*: foliar magnesium (g<sup>1</sup> g<sup>-1</sup>)
- *Al_foliar*: foliar aluminum (g<sup>1</sup> g<sup>-1</sup>)
- *Zn_foliar*: foliar zinc (g<sup>1</sup> g<sup>-1</sup>)
- *Ca_soil*: soil calcium (g<sup>1</sup> g<sup>-1</sup>)
- *P_soil*: soil phosphorus (g<sup>1</sup> g<sup>-1</sup>)
- *K_soil*: soil potassium (g<sup>1</sup> g<sup>-1</sup>)
- *Mg_soil*: soil magnesium (g<sup>1</sup> g<sup>-1</sup>)
- *Al_soil*: soil aluminum (g<sup>1</sup> g<sup>-1</sup>)
- *Zn_soil*: soil zinc (g<sup>1</sup> g<sup>-1</sup>)
- *pH*: soil pH
- *CEC*: soil cation exchange capacity (cmol<sub>c</sub> kg<sup>-1</sup>)
- *C_soil*: soil carbon (g<sup>1</sup> g<sup>-1</sup>)
- *N_soil*: soil nitrogen (g<sup>1</sup> g<sup>-1</sup>)
- *Retention*: soil water retention (%)

### /data/mdi_stand_density.csv
This file contains stand density data. The specific column information is:
- 

### /scripts/mdi_pitchpine_analyses.csv
This file contains R code with the data analysis, graph and table creation.
