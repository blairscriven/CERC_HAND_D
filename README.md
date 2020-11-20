# CERC_HAND_D
Canadian Estimator of Ratings Curves using HAND and Discharge (CERC-HAND-D). A custom ArcGIS Pro tool used to create a synthetic rating curve using publicly available data from the  Government of Canada.

This tool was designed to work with input data soured from the Open Government database provided by the Government of Canada. The suggested data to use includes:

2015 Land Cover of Canada: https://open.canada.ca/data/en/dataset/4e615eae-b90c-420b-adee-2ca35896caf6

High Resolution Digital Elevation Model (HRDEM) - CanElevation Series: https://open.canada.ca/data/en/dataset/957782bf-847c-4644-a757-e383c0057995

National Hydrographic Network: https://www.nrcan.gc.ca/science-and-data/science-and-research/earth-sciences/geography/topographic-information/geobase-surface-water-program-geeau/national-hydrographic-network/21361

CERC-HAND-D requires a Height Above Nearest Drainage (HAND) raster file; a HAND raster can be created using the datasets listed above. Alternative data sources may reduce the quality of the results or prevent the tool from functioning. The National HydroGraphic Network does not provide attribute data on Stahler order, so an alternative source for river network data may be preferred.Ideally, the area of interest should have around 1.5-3 km length of river (high order stahler) and a river 
slope less than 0.001%. Keeping within these guidelines will produce the most accurate results possible. 
