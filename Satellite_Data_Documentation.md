## Satellite Data for FRESCA Project

Satellite-based datasets for the FRESCA project:
USF has compiled several satellite-based datasets, including ocean color, sea surface temperature and biogeographic seascapes for waters around South Florida. A Github repository was created to share code among the project team. Details on data products along with links to data and the Github repository are given below.

#### Data Types (details on each are given below):
1. Ocean color: from MODIS-Aqua, VIIRS-SNPP and PACE-OCI sensors
2. Sea Surface Temperature (SST): from Multi-scale Ultra-high Resolution (MUR) SST Analysis fv04.1 and VIIRS-SNPP
3. Biogeographic Seascapes: Based on MODIS-Aqua data and served by NOAA Coastwatch
4. Net Primary Productivity: From Oregon State University Ocean Productivity Lab

#### NOTE on MODIS-Aqua: MODIS-Aqua is at the end of its operational life and has exhibited degraded data quality since the end of 2022. While there have been data reprocessings bu OB-DAAC to improve data quality, MODIS will be decommissioned at some point in 2026. We will continue to serve these products with reprocessed data as long as they are available due to the long MODIS data record. Users should plan to transition to using products from VIIRS-SNPP and PACE-OCI. Please feel free to contact Dan Otis at USF with questions: dotis@usf.edu

ERDDAP Link for all datasets:  
http://131.247.136.200:8080/erddap/info/index.html

### OCEAN COLOR
Data are obtained from NASA's Ocean Biology Distributed Active Archive Center (OB-DAAC):  
https://oceancolor.gsfc.nasa.gov/

Data is ingested, composited to daily, weekly and monthly time intervals, gridded and served via ERDDAP at USF.

Within each Ocean Color file, there are the following products:
1. chlor_a (chlorophyll-a concentration in mg m^-3) - proxy for phytoplankton biomass
2. chlor_a_clim (climatology of chlorophyll-a concentration based on data from 2003-2019 for MODIS-Aqua and 2013-2019 for VIIRS)
3. chlor_a_anom (anomaly of chlorophyll-a concentration; chlor_a_clim subtracted from chlor_a)
4. Rrs_667 (remote sensing reflectance at 667nm in sr ^-1) - proxy for suspended sediments
5. Rrs_667_clim (climatology of Rrs_667 based on data from 2003-2019 for MODIS-Aqua and 2013-2019 for VIIRS)
6. Rrs_667_anom (anomaly of Rrs_667; Rrs_667_clim subtracted from Rrs_667)
7. Kd_490 (diffuse attenuation coefficient at 490nm in m^-1) - proxy for water clarity
8. Kd_490_clim (climatology of Kd_490 based on data from 2003-2019 for MODIS-Aqua and 2013-2019 for VIIRS)
9. Kd_490_anom (anomaly of Kd_490; Kd_490_clim subtracted from Kd_490)
##### The following products are included with MODIS-Aqua files only
10. ABI (Algal Bloom Index in mW cm^2 um^-1 sr^-1) - proxy for potential algal bloom conditions
11. ABI_clim (climatology of ABI based on data from 2003-2019 for MODIS-Aqua and 2013-2019 for VIIRS)
12. ABI_anom (anomaly of ABI; ABI_clim subtracted from ABI)

#### List of Ocean Color Data Products
1. MODIS-Aqua Ocean Color 7-Day means for GOM (18N to 31N; -98E to -78.5E)  
ERDDAP Link: http://131.247.136.200:8080/erddap/griddap/moda_oc_7d_gom.graph

2. MODIS-Aqua Ocean Color Monthly means for GOM  (18N to 31N; -98E to -78.5E)  
ERDDAP Link: http://131.247.136.200:8080/erddap/griddap/moda_oc_mo_gom.graph

3. MODIS-Aqua Ocean Color 1-Day composite for Florida (24N to 31N; -85E to -78.5E)  
ERDDAP Link: MODIS-AQUA Ocean Color 1-Day composite for Florida

4. MODIS-Aqua Ocean Color 7-Day means for the Southeast US (29N to 40.5N; -82E to -73E)  
ERDDAP Link: http://131.247.136.200:8080/erddap/griddap/moda_oc_7d_seus.graph

5. VIIRS-SNPP Ocean Color 7-Day means for GOM (18N to 31N; -98E to -78.5E)  
ERDDAP Link: http://131.247.136.200:8080/erddap/griddap/vsnpp_oc_7d_gom.graph

6. VIIRS-SNPP Ocean Color Monthly means for GOM  (18N to 31N; -98E to -78.5E)  
ERDDAP Link: http://131.247.136.200:8080/erddap/griddap/vsnpp_oc_mo_gom.graph

7. VIIRS-SNPP Ocean Color 1-Day composite for Florida (24N to 31N; -85E to -78.5E)  
ERDDAP Link: http://131.247.136.200:8080/erddap/griddap/vsnpp_oc_1d_fl.graph

8. VIIRS-SNPP Ocean Color 1-Day composite for NW Gulf of Mexico (24N to 31N; -85E to -78.5E)  
ERDDAP Link: http://131.247.136.200:8080/erddap/griddap/vsnpp_oc_1d_nwgom.graph

9. VIIRS-SNPP Ocean Color 7-Day means for the Southeast US (29N to 40.5N; -82E to -73E)  
ERDDAP Link: http://131.247.136.200:8080/erddap/griddap/vsnpp_oc_7d_seus.graph

10. PACE-OCI OC AOP (incl. hyperspectral Rrs and nflh) 1-Day composite for South Florida (24N to 31N; -85E to -78.5E)
ERDDAP Link: http://131.247.136.200:8080/erddap/griddap/pace_oc_aop_1d_sfl.graph
 - These files contain remote sensing reflectance at all visible PACE-OCI bands and normalized fluorescence line height

11. PACE-OCI Ocean Color BGC 1-Day composite for South Florida (24N to 31N; -85E to -78.5E)
ERDDAP Link: http://131.247.136.200:8080/erddap/griddap/pace_oc_bgc_1d_sfl.graph
 - These files contain chlorophyll-a concentration and phytoplankton carbon from PACE-OCI

### SEA SURFACE TEMPERAURE
Sea surface temperature products come from two main sources:
1. Multi-scale Ultra-high Resolution (MUR) SST Analysis fv04.1 (https://podaac.jpl.nasa.gov/dataset/MUR-JPL-L4-GLOB-v4.1)
 - This is a global dataset produced as part of the Global High Resolution Sea Surface Temperature (GHRSST) project
 - It is based on data from a large suite of sensors
2. VIIRS-SNPP 1-km nighttime SST
 - These products are obtained from the OB-DAAC and composited, gridded and served via ERDDAP similar to the Ocean Color products listed above

#### List of SST Data Products
1. 	Multi-scale Ultra-high Resolution (MUR) SST Analysis fv04.1, Global, 0.01°, 2002-present, Daily (-90N to 90N; -180E to 180E)
ERDDAP Link: http://131.247.136.200:8080/erddap/griddap/jplMURSST41.graph

2. Multi-scale Ultra-high Resolution (MUR) SST Analysis Anomaly fv04.1, Global, 0.01°, 2002-present, Daily (-90N to 90N; -180E to 180E)  
ERDDAP Link: http://131.247.136.200:8080/erddap/griddap/jplMURSST41anom1day.graph

3. VIIRS-SNPP Nighttime SST 7-Day mean for GOM  (24N to 31N; -85E to -78.5E)    
ERDDAP Link: http://131.247.136.200:8080/erddap/griddap/vsnpp_sstn_7d_gom.graph

4. 	VIIRS-SNPP Nighttime SST 7-Day mean for SEUS (29N to 40.5N; -82E to -73E)  
ERDDAP Link: http://131.247.136.200:8080/erddap/griddap/vsnpp_sstn_7d_seus.graph


### DYNAMIC BIOGEOGRAPHIC SEASCAPES
Dynamic seascapes are classified using multivariate satellite data and describe basin and gyre scale features and seasonal boundary shifts. They are obtained from [NOAA's Coastwatch](https://cwcgom.aoml.noaa.gov/erddap/griddap/noaa_aoml_seascapes_8day.graph). A description of seascape classes can be found [HERE](https://shiny.marinebon.app/seascapes/classes.html).
The seascapes served here are based on MODIS-Aqua data (-90N to 90N; -180E 180E)  
ERDDAP Link: http://131.247.136.200:8080/erddap/griddap/noaa_aoml_seascapes_8day.graph

There is a version based on VIIRS-SNPP and Sentinel-3 OLCI data that is served by NOAA Coastwatch  
ERDDAP Link: https://cwcgom.aoml.noaa.gov/erddap/griddap/noaa_aoml_seascapes_8day_viirs_olci.graph


### NET PRIMARY PRODUCTIVITY
Net primary productivity products are calculated using the Vertically Generalized Productivity Model (VGPM) based on MODIS-Aqua data at 9-km spatial resolution (-90N to 90N; -180E to 180E)  
These products are from the Oregon State Ocean Productivity Lab (http://orca.science.oregonstate.edu/index.php)  
ERDDAP Link: http://131.247.136.200:8080/erddap/griddap/moda_npp_mo_glob.graph

