### Creates a Synthetic Ratings Table (relationship between discharge (Q) and
### stage height(H)) for Canadian rivers. Uses river geometry
###  and Manning's equation to find the relationship; two columns are created,
### one using the minimum roughness coefficient value, and the other using 
### the median value (the maximum roughness coefficient value was experimented
### with, but ultimately was found to highly overestimate flood levels). The 
### minimum and median values create a minimum and maximum flood scenarios,
### which in turn creates prediction window. 
### Requires a Height Above nearest Drainage (HAND) raster, a single polygon 
### shapefile that encompasses the entire watershed, the 2015 Land Cover of 
### Canada file, a DEM or DTM file, and stahler high order stream network.
### The code is based off of equations 5-11
### from Zheng, X., Tarboton, D. G., Maidment, D. R., Liu, Y. Y.,
### & Passalacqua, P. (2018). River channel geometry and rating
### curve estimation using height above the nearest drainage. JAWRA
### Journal of the American Water Resources Association, 54(4), 785-806.

import arcpy
from arcpy import *
import arcpy.sa as SA
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Set the zonal shape file needed for all Zonal Statistics
inZoneData = arcpy.GetParameterAsText(2)
zoneField = "FID"

#Insert Hand Model and Stage  level
HAND_RAW = arcpy.GetParameterAsText(0)
HAND = SA.ExtractByMask(HAND_RAW, inZoneData)
MAXwaterLevel = float(arcpy.GetParameterAsText(1))

#Define variables needed to run the flood simumlations loops
Water_Level = 0.1 # Base Water lvl iteration
loop = 0

#Input River length
All_River_lines = arcpy.GetParameterAsText(3)
Data_output_Location = arcpy.GetParameterAsText(4) #This will be used to create interprocessing shapefiles

Output_Flowlines = Data_output_Location + "/Flowlines.shp"
arcpy.Intersect_analysis ([inZoneData, All_River_lines], Output_Flowlines, "ONLY_FID", "", "line")

arcpy.AddField_management(Output_Flowlines,'LENGTH','DOUBLE')
arcpy.CalculateField_management(Output_Flowlines,'LENGTH','!shape.length@meters!','PYTHON')



RiverLength = int(0)

with arcpy.da.SearchCursor(Output_Flowlines, 'SHAPE@LENGTH') as cursor:
    for row in cursor:
        RiverLength = RiverLength + row[0]


### Start OF Get Min-Medium Roughness CODE
### Racll that we use medium instead of max (from experiments)

### Insert the LandClass, run the Max flood level on the HAND Model, Cut out 
### the overlapping LandClass cells, then add all the Manning's Coefficient
LandClass_Raw = arcpy.GetParameterAsText(5)
MaxHand = float(MAXwaterLevel) - Raster(HAND)
MaxHand_rast = SA.Con(MaxHand, MaxHand, "", "VALUE > 0")
LandClass = SA.ExtractByMask(LandClass_Raw, MaxHand_rast)

#create Landclass table so that cursor (update and search) can be used
arcpy.TableToTable_conversion(LandClass, Data_output_Location, "LandClass_Table.dbf")
LandClass_table = os.path.join(Data_output_Location, "LandClass_Table.dbf")

#Set the roughness for the river channel
Water_Rough_text = arcpy.GetParameterAsText(7)
Water_Rough = 0.045
if Water_Rough_text == "Clean; Straight":
    Water_Rough = 0.03
elif Water_Rough_text == "Rocks+Vegetation; Straight":
    Water_Rough = 0.035
elif Water_Rough_text == "Clean; Winding":
    Water_Rough = 0.04
elif Water_Rough_text == "Rocks+Vegetation; Winding":
    Water_Rough = 0.045
Water_Rough_flt = float(Water_Rough)
arcpy.AddMessage("Channel Rough: " + str(Water_Rough))

arcpy.AddField_management(LandClass_table,'RoughCoef','DOUBLE')
Roughness_update_Feilds = ['Value', 'RoughCoef']
with arcpy.da.UpdateCursor(LandClass_table, Roughness_update_Feilds) as cursor_area:
    for row in cursor_area:
       options_dict = {1: 0.1, 2: 0.1, 5: 0.1, 6: 0.1, 8: 0.07, 10: 0.03, 11: 0.05,
                         12: 0.035, 13: 0.027, 14: 0.035, 15: 0.037, 16: 0.03, 17: 0.05,
                          18: Water_Rough, 19: 0.07}
       row[1] = options_dict.get(row[0], 0.033)
       cursor_area.updateRow(row)



### Get Min(Q1) and Medium (Q2) roughness values
Landfields_Q1 = ['RoughCoef']
roughlist = []
rows = arcpy.da.SearchCursor(LandClass_table, Landfields_Q1)
for row in rows:
    val = float(row[0])
    roughlist.append(val)

Rough_Q1 = min(roughlist)

roughlist.sort() 
mid = len(roughlist) // 2
Rough_Q2 = (roughlist[mid] + roughlist[~mid]) / 2

#Convert all necessary variables into float for the final discharge calculation
Rough_Q1_flt = float(Rough_Q1)
Rough_Q2_flt = float(Rough_Q2)

### END OF Get Min-Medium Roughness CODE


###Now to get all the other Manning equation inputs and then calculate Stream Discharge
# Insert Digital elevation Model (raw)
DEM_RAW = Raster(arcpy.GetParameterAsText(6))
DEM = SA.ExtractByMask(DEM_RAW, inZoneData)

###### CODE TO FIND RIVER SLOPE
# Create dangling points, using Feature Vertices To Points tool, on the clipped 
# section of the hydro-network (within the watershed)
Output_Verticles = Data_output_Location + "/verticles.shp"
arcpy.FeatureVerticesToPoints_management(Output_Flowlines,
                                         Output_Verticles, 
                                         "DANGLE")

# Find elevation values at dangling points using Extract Values to Points tool.
# Creates column called "RASTERVALU" that contains elevation samples
Output_Verticles_elevation = Data_output_Location + "/verticles_Elevation.shp"
arcpy.sa.ExtractValuesToPoints(Output_Verticles, DEM_RAW,
                      Output_Verticles_elevation,"INTERPOLATE",
                      "VALUE_ONLY")

# Now to put the elevation samples into a list (Elevation_list)
Elevation_column = ['RASTERVALU']
Elevation_list = []
Elevation_rows = arcpy.da.SearchCursor(Output_Verticles_elevation, Elevation_column)
for Elevation_row in Elevation_rows:
    Elevation_val = float(Elevation_row[0])
    Elevation_list.append(Elevation_val)

# Find the maximum and minimum elevation samples
Maximum_Elevation = max(Elevation_list)
Minimum_Elevation = min(Elevation_list)

Slope = (Maximum_Elevation - Minimum_Elevation)/ RiverLength
if Slope > 1 or Slope <= 0:
    Slope = 0.002



arcpy.AddMessage("Slope: " + str(Slope))

arcpy.Delete_management(Output_Verticles)
arcpy.Delete_management(Output_Verticles_elevation)
###### END OF CODE TO FIND RIVER SLOPE

#Overland Water Level Simulation iterations
iteration = float(arcpy.GetParameterAsText(8))

# Execute CreateTable (Create the Flood-Discharge Table)
tableName = "RatingCurve_Manning_Data.dbf"
arcpy.CreateTable_management(Data_output_Location, tableName)

# Set variables for AddField
out_name = os.path.join(Data_output_Location, tableName)
fieldName1 = "WaterLvl_H"
fieldName2 = "ChanWidth"
fieldName3 = "CrossArea"
fieldName4 = "WetPerm"
fieldName5 = "HydroRad"
fieldName6 = "RivLength"
fieldName7 = "Q_Min"
fieldName8 = "Q_Max"

# Execute AddField 8 times for 8 new fields
arcpy.AddField_management(out_name, fieldName1, "DOUBLE", 10, 7, field_is_nullable="NON_NULLABLE")
arcpy.AddField_management(out_name, fieldName2, "DOUBLE", 10, 7, field_is_nullable="NON_NULLABLE")
arcpy.AddField_management(out_name, fieldName3, "DOUBLE", 10, 7, field_is_nullable="NON_NULLABLE")
arcpy.AddField_management(out_name, fieldName4, "DOUBLE", 10, 7, field_is_nullable="NON_NULLABLE")
arcpy.AddField_management(out_name, fieldName5, "DOUBLE", 10, 7, field_is_nullable="NON_NULLABLE")
arcpy.AddField_management(out_name, fieldName6, "DOUBLE", 10, 7, field_is_nullable="NON_NULLABLE")
arcpy.AddField_management(out_name, fieldName7, "DOUBLE", 10, 7, field_is_nullable="NON_NULLABLE")
arcpy.AddField_management(out_name, fieldName8, "DOUBLE", 10, 7, field_is_nullable="NON_NULLABLE")


# Execute DeleteField to remove default field created when CreateTable was executed.
arcpy.DeleteField_management(out_name, ["Field1"])




while Water_Level <= MAXwaterLevel:
    
    loop = loop + 1    
    arcpy.AddMessage("Simulating Overflow Water Level: " + str(Water_Level) + "m")
    #Calculate Depths
    depth = float(Water_Level) - Raster(HAND)
    
    #Conditional statement to replace negative values or 0 with nodata; these values do not show flood extent
    rast = SA.Con(depth, depth, "", "VALUE > 0") 


    #Calculating cell area; not the total area
    descRast = arcpy.Describe(rast)
    x_cell = descRast.meanCellWidth
    y_cell = descRast.meanCellHeight
    cellArea = x_cell * y_cell


    #Calculate pixel count and Surface Area using Zonal Statistics as Table
    #and a cursor to extract those values
    outTable = Data_output_Location + "/ZoneStatsTable.dbf"
    outZSaT = SA.ZonalStatisticsAsTable(inZoneData, zoneField, rast, outTable, "DATA", "SUM")

    fc = Data_output_Location + "/ZoneStatsTable.dbf"
    field_count = "COUNT"
    feild_area = "AREA"
    cursor_area = arcpy.SearchCursor(fc)
    for row in cursor_area:
        pixelCount = int(row.getValue(field_count))
        SurfaceArea =int(row.getValue(feild_area))



    #Calculating Flood Volume of inundation zone; eqn 6 from Zheng et al.(2018)
    FloodVolRast = rast * cellArea 
    FloodVolZonalStatistics = SA.ZonalStatistics(inZoneData, zoneField, FloodVolRast, "SUM", "DATA")
    FloodVolMaxResult = arcpy.GetRasterProperties_management(FloodVolZonalStatistics, "MAXIMUM")
    FloodVol = float(FloodVolMaxResult.getOutput(0))
    FloodVol_str = str(FloodVol)
    print("Flood Volune is (and it should equal test_table): " + FloodVol_str + "m^3")



    #Calculating channel bed area of inundation zone; eqn 5 from Zheng et al. (2018)
    DEM_cut = (DEM + rast) - rast
    DEMslope = SA.Slope(DEM_cut, "PERCENT_RISE")
    DEMRiseRun = DEMslope / 100 
    ChannelBedRast = cellArea * SA.SquareRoot(1 + SA.Square(DEMRiseRun))
    ChannelBedZonalStatistics = SA.ZonalStatistics(inZoneData, zoneField, ChannelBedRast, "SUM", "DATA")
    ChannelBedMaxResult = arcpy.GetRasterProperties_management(ChannelBedZonalStatistics, "MAXIMUM")
    ChannelBed = float(ChannelBedMaxResult.getOutput(0))
    ChannelBed_str = str(ChannelBed)
    print("Channel bed area: " + ChannelBed_str + "m^2")
    


    ###Calculate for eqn 7-10 in Zheng et al. (2018)

    #Reach-average channel width; eqn 7 from Zheng et al. (2018)
    Channel_width = SurfaceArea / RiverLength
    Channel_width_str = str(Channel_width)
    print("Reach-average channel width = " + Channel_width_str + " m")
    Channel_width_flt = float(Channel_width)

    #Reach-average cross section area; eqn 8 from Zheng et al. (2018)
    cross_section_area = FloodVol / RiverLength # Could add river base (height * width of river at normal) 
    cross_section_area_str = str(cross_section_area)
    print("Reach-average cross section area = " + cross_section_area_str + " m^2")
    cross_section_area_flt = float(cross_section_area)

    #Reach-average cross section wetted perimeter; eqn 9 from Zheng et al. (2018)
    cross_section_perm = ChannelBed / RiverLength
    cross_section_perm_str = str(cross_section_perm)
    print("Reach-average cross section wetted perimeter = " + cross_section_perm_str + " m")
    cross_section_perm_flt = float(cross_section_perm)

    #Hydraulic Radius; eqn 10 from Zheng et al. (2018)
    HydroRadius = cross_section_area / cross_section_perm
    HydroRadius_str = str(HydroRadius)
    print("Hydraulic Radius = " + HydroRadius_str + " m")
    
    #Area_flt = float(Area)
    Slope_flt = float(Slope)
    HydroRadius_flt = float(HydroRadius)
    
    ###Now to calculate all three of the Discharges
    #Flood Discharge equation (Rough_Q1)
    Q_min = (1.0 / Rough_Q1) * (HydroRadius_flt ** 0.666) * cross_section_area_flt * (Slope_flt ** 0.5)

    #Flood Discharge equation (Rough_Q2)
    Q_med = (1.0 / Rough_Q2) * (HydroRadius_flt ** 0.666) * cross_section_area_flt * (Slope_flt ** 0.5)

    #Now to convert all the remaining values to float so that they can be added to the FloodDischarge table
    Q_min_flt = float(Q_min)
    Q_med_flt = float(Q_med)
    Rivlength_flt = float(RiverLength) #This is here because I didn't know where else to put it



     ###Code to add values to FloodDischarge table


    # Set variables for InsertCursor
    fieldNames = [fieldName1, fieldName2, fieldName3, fieldName4, fieldName5, fieldName6, fieldName7, fieldName8]


    row_values = [(Water_Level, Channel_width_flt, cross_section_area_flt, cross_section_perm_flt,
                                HydroRadius_flt, Rivlength_flt, Q_min_flt, Q_med_flt)]

    # Open an InsertCursor
    cursor_discharge = arcpy.da.InsertCursor(out_name, fieldNames)

    # Insert new rows
    for row in row_values:
        cursor_discharge.insertRow(row)

    # Delete cursor
    del cursor_discharge
    
    # Keep the loop going
    Push_lvl = iteration * (loop)
    Water_Level = round(Push_lvl, 2)

arcpy.Delete_management(outTable)
arcpy.Delete_management(LandClass_table)

###CODE to create the log and power figures!
arcpy.AddMessage("Creating Figures....")
Qmin = [] # Min n Discharge
Qmed = [] # Med n Discharge
H = [] # Stage Heights


with arcpy.da.SearchCursor(out_name, fieldNames) as cursor:
    for row in cursor:
        H.append(row[0])
        Qmin.append(row[6])
        Qmed.append(row[7])

xmin = np.array(Qmin) # Min n Discharge
xmed = np.array(Qmed) # Med n Discharge
y = np.array(H) # Stage Heights

# find the variables in the log and power equations
logarithmic_equation_min = np.polyfit(np.log(xmin), y, 1)
power_equation_min, pcov = curve_fit(lambda fx,a,b: a*fx**-b,  xmin,  y)
logarithmic_equation_med = np.polyfit(np.log(xmed), y, 1)
power_equation_med, pcov = curve_fit(lambda fx,a,b: a*fx**-b,  xmed,  y)

# Define the log and power equations
ylog_min = (logarithmic_equation_min[0] * np.log(xmin)) + logarithmic_equation_min[1]
ypwr_min = power_equation_min[0]*xmin**-power_equation_min[1]
ylog_med = (logarithmic_equation_med[0] * np.log(xmed)) + logarithmic_equation_med[1]
ypwr_med = power_equation_med[0]*xmed**-power_equation_med[1]

# Create Labels for the graphs
log_min_legendlabel = "H = " + str(round(logarithmic_equation_min[0], 3)) + "ln(Q) + (" + str(round(logarithmic_equation_min[1], 3)) + ")"
pwr_min_legendlabel = "H = " + str(round(power_equation_min[0], 3)) + " * Q^-(" + str(round(power_equation_min[1], 3)) + ")"
log_med_legendlabel = "H = " + str(round(logarithmic_equation_med[0], 3)) + "ln(Q) + (" + str(round(logarithmic_equation_med[1], 3)) + ")"
pwr_med_legendlabel = "H = " + str(round(power_equation_med[0], 3)) + " * Q^-(" + str(round(power_equation_med[1], 3)) + ")"


# Code to create the graphs
graph, node = plt.subplots(1,2, figsize=(15,5))
node[0].set_title("Power Curve - Prediction Window")
node[0].scatter(xmin, y, c='green')
node[0].plot(xmin, ypwr_min, c= "green", linestyle='solid', label = pwr_min_legendlabel)
node[0].scatter(xmed, y, c='red')
node[0].plot(xmed, ypwr_med, c= "red", linestyle='solid', label = pwr_med_legendlabel)
node[0].set_xlabel("Discharge / Q (m^3/sec)")
node[0].set_ylabel("Water Level / H (m)")
node[0].legend(loc="upper left")

node[1].set_title("Logarithmic Curve - Prediction Window")
node[1].scatter(xmin, y, c='green')
node[1].plot(xmin, ylog_min, c= "green", linestyle='solid', label = log_min_legendlabel)
node[1].scatter(xmed, y, c='red')
node[1].plot(xmed, ylog_med, c= "red", linestyle='solid', label = log_med_legendlabel)
node[1].set_xlabel("Discharge / Q (m^3/sec)")
node[1].set_ylabel("Water Level / H (m)")
node[1].legend(loc="lower right")


plt.show()
###END OF CODE to create the log and power figures!
