# Contains script and some data tables for FRESCA Module 1 multi-stressor analysis

## **Meeting notes:**

**May 14th 2024: Dan Otis and Enrique Montes**
- Daily and 7-day 1 km MODIS and VIIRS data is available on the [project ERDDAP server](http://131.247.136.200:8080/erddap/info/index.html?page=1&itemsPerPage=1000) for:
  - Aqua chlorophyll-a
  - nFLH
  - turbidity (Rrs 667)
  - Alagae Bloom Index (nFLH/Rrs547)
  - SST
  - SST anomaly 
- Seagrass maps:
  - we should use what we already have at this point: map with presence/absense and maybe a time-series plot with seagrass area for Florida Bay.
  - Carolina Peralta will take over this portion of the project
  - define priority areas looking forward to generate new maps and time series analysis
  - Carolina should generate these products. This requires:
    - Selecting LandSat and Sentinel RGBs for quality control by hand
    - Validation with in situ data, which we have for FL Bay
    - Best images around the winter-spring time (February through April). Summertime is not good due to high humidity and turbidity from runoff.
  - Frank, Dan, Carolina and Enrique will meet after first-year reporting to define next steps.  



**February 9th 2024: Dan Otis and Enrique Montes**
- We will use [Multi-scale Ultra-high Resolution (MUR) SST Analysis](https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.html) for thermal stress distributions
- Dan will upload 1 km MODIS Aqua chlorophyll-a, nFLH, turbidity, and Alagae Bloom Index data for the entire Gulf of Mexico (2002-present) to an ERDDAP server
- Dan will compile freshwater discharge data for: Caloosahatchee River, Shark River, Manatee River, Peace River, and maybe other.
- Goals:
  - Spatially characterize stressors
  - Generate products showing extent and timing of stressor events: heatwaves, cold spells, freshwater discharge, red tides, others
  - Define stressor thresholds that can be used for model parameterization
  - Match ABI data with red tide records
