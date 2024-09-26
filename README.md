# Grazing-niche

Code to reproduce the analysis and figures for the article:
CH. Li, M. Kotz, P. Pradhan, XD. Wu, Z. Li, GQ. Chen. Global Climate Niche for Grazing Implies Shifting Suitability
 
 
1. Data

Primary data are not provided in the repository due to storage capacity. Full download of these data requires storage capacity of 246GB. All primary data used in this study are openly accessible at the following sites:
 

CMIP (Coupled Model Intercomparison Project) can be accessed at https://pcmdi.llnl.gov/CMIP5/. 
 

ESFG data (Earth System Grid Federation) can be accessed at https://esgf.llnl.gov CHELSA V2 (Climatologies at high resolution for the earth’s land surface areas) can be accessed at https://chelsa-climate.org. 
 

Hyde (History database of the Global Environment) can be accessed at https://www.pbl.nl/en/hyde-history-database-of-the-global-environment. 
 

Gridded Livestock of the World (GLW 3) database can be accessed at https://livestockdata.org/contributors/food-and-agriculture-organization-united-nations/gridded-livestock-world-glw3. 
 

Gridded Population of the World Version 4 (GPWv4) from SEDAC (Socioeconomic Data and Applications Center) can be accessed at https://sedac.ciesin.columbia.edu/data/collection/gpw-v4. 
 

The MOD17A3HGF Version 6 dataset can be accessed at https://lpdaac.usgs.gov/products/mod17a3hgfv006/. 
 

Percent Tree Coverage (PTC) Global version can be accessed at https://globalmaps.github.io/ptc.html. 
 

EARTH Env data can be accessed at https://www.earthenv.org. Land-Use Harmonization (LUH2) can be accessed at https://luh.umd.edu. 
 

Global Aboveground and Belowground Biomass Carbon Density Maps can be accessed at https://daac.ornl.gov/VEGETATION/guides/Global_Maps_C_Density_2010.html. 
 

SOC data is from literature (https://pnas.org/doi/full/10.1073/pnas.1706103114. 
 

2. Code for replication
 

Code for historical climate data for GISS, historical pasture data, niche threshold for global and continental values are in script code_climate1_land. Analysis for two new climate models are in code_climate2.m. Analysis for 4 extra climate models are run at code_climate_data.m
 

Code to run present climate data and land use, rescaling to align with historical data are in  script code_nicheanalysis.m
 

Code to run analysis for livestock distribution, above-ground biomass carbon, below ground biomass carbon, SOC, are in script code_livestock_distribution
 

Code to process data from future land use, future climate change under rcp scenarios, future niche threshold, are run at code_future_climate.m
 

Code to run overgrazing analysis are run in code_overgraze_analysis.m


Code to run uncertainty analysis for livestock distribution, carbon storage are in script code_montecarlo.m
 

Code to run analysis for future change of grassland, and future climate change, future shift of niche range, and future continental analysis are run at code_future_grassland.m 
 

Code to run reanalysis by using smaller niche range are in script code_regional_analysis.m and also code_climate1_land
 

Code for analyzing niche for regional scale are in script code_regional_climate.m
 

Code to reproduce the figures are run at code_figures.m
 

3. Dependencies required
 

For Matlab, our scripts require the dependencies require the climate toolbox
 
(https://zenodo.org/badge/latestdoi/171331090), land mask, and addcolorplus package.
 

All scripts were run in version Matlab_R2021b environment.
 

For any further questions please feel free to contact: Li.Chaohui@pik-potsdam.de
 
