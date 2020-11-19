### Creates a Synthetic Ratings Table (relationship between discharge and
### stage height) for Canadian rivers (High order stahler). Uses river geometry
###  and Manning's equation to find the relationship; three columns are created,
### one using the minimum roughness coefficient value, and the other two using 
### the median and maximum value. 
### Requires a Height Above nearest Drainage (HAND) raster, a single polygon 
### shapefile that encompasses the entire watershed, the 2015 Land Cover of 
### Canada file, a DEM or DTM file, and high order stream network.
### The code is based off of equations 5-11
### from Zheng, X., Tarboton, D. G., Maidment, D. R., Liu, Y. Y.,
### & Passalacqua, P. (2018). River channel geometry and rating
### curve estimation using height above the nearest drainage. JAWRA
### Journal of the American Water Resources Association, 54(4), 785-806.

import arcpy
from arcpy import *
import arcpy.sa as SA
import os

# Set the zonal shape file needed for all Zonal Statistics
inZoneData = arcpy.GetParameterAsText(4)
zoneField = "FID"

#Insert Hand Model and Stage  level
HAND_RAW = arcpy.GetParameterAsText(0)
HAND = SA.ExtractByMask(HAND_RAW, inZoneData)
MAXwaterLevel = float(arcpy.GetParameterAsText(1))

#Define water level iteration for the LOOP
Base_River_Length = float(arcpy.GetParameterAsText(2))
Base_River_Depth = float(arcpy.GetParameterAsText(3))
Base_River_Volume = (Base_River_Length * Base_River_Depth) / 2
Water_Level = 0


#Input River length
All_River_lines = arcpy.GetParameterAsText(5)
Data_output_Location = arcpy.GetParameterAsText(6) #This will be used to create interprocessing shapefiles

Output_Flowlines = Data_output_Location + "/Flowlines.shp"
arcpy.Intersect_analysis ([inZoneData, All_River_lines], Output_Flowlines, "ONLY_FID", "", "line")

arcpy.AddField_management(Output_Flowlines,'LENGTH','DOUBLE')
arcpy.CalculateField_management(Output_Flowlines,'LENGTH','!shape.length@meters!','PYTHON')



RiverLength = int(0)

with arcpy.da.SearchCursor(Output_Flowlines, 'SHAPE@LENGTH') as cursor:
    for row in cursor:
        RiverLength = RiverLength + row[0]


### Insert the LandClass, run the Max flood level on the HAND Model, Cut out 
### the overlapping LandClass cells, then add all the Manning's Coefficient

LandClass_Raw = arcpy.GetParameterAsText(7)
MaxHand = float(MAXwaterLevel) - Raster(HAND)
MaxHand_rast = SA.Con(MaxHand, MaxHand, "", "VALUE > 0")
LandClass = SA.ExtractByMask(LandClass_Raw, MaxHand_rast)

#create Landclass table so that cursor (update and search) can be used
arcpy.TableToTable_conversion(LandClass, Data_output_Location, "LandClass_Table.dbf")
LandClass_table = os.path.join(Data_output_Location, "LandClass_Table.dbf")

#Set the roughness for the river channel
Water_Rough_text = arcpy.GetParameterAsText(9)
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
        if row[0] == 1 :
            row[1] = 0.1
        elif row[0] == 2 :
            row[1] = 0.1
        elif row[0] == 5 :
            row[1] = 0.1
        elif row[0] == 6 :
            row[1] = 0.1
        elif row[0] == 8 :
            row[1] = 0.07
        elif row[0] == 10 :
            row[1] = 0.03
        elif row[0] == 11 :
            row[1] = 0.05
        elif row[0] == 12 :
            row[1] = 0.035
        elif row[0] == 13 :
            row[1] = 0.027
        elif row[0] == 14 :
            row[1] = 0.035
        elif row[0] == 15 :
            row[1] = 0.037
        elif row[0] == 16 :
            row[1] = 0.03
        elif row[0] == 17 :
            row[1] = 0.05
        elif row[0] == 18 :
            row[1] = Water_Rough
        elif row[0] == 19 :
            row[1] = 0.07
        else:
            row[1] = 0.033
        cursor_area.updateRow(row)




### Get Min-Medium-Max
Landfields_Q1 = ['RoughCoef']
roughlist = []
rows = arcpy.da.SearchCursor(LandClass_table, Landfields_Q1)
for row in rows:
    val = float(row[0])
    roughlist.append(val)


###Now to get all the other Manning equation inputs and then calculate Stream Discharge
Rough_Q1 = min(roughlist)
Rough_Q3 = max(roughlist)

roughlist.sort() 
mid = len(roughlist) // 2
Rough_Q2 = (roughlist[mid] + roughlist[~mid]) / 2

#Convert all necessary variables into float for the final discharge calculation
Rough_Q1_flt = float(Rough_Q1)
Rough_Q2_flt = float(Rough_Q2)
Rough_Q3_flt = float(Rough_Q3)


### END OF Get Min-Medium-Max CODE



# Insert Digital elevation Model (raw)
DEM_RAW = Raster(arcpy.GetParameterAsText(8))
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

arcpy.AddMessage("Maximum_Elevation: " + str(Maximum_Elevation))
arcpy.AddMessage("Minimum_Elevation: " + str(Minimum_Elevation))
arcpy.AddMessage("Slope: " + str(Slope))
arcpy.AddMessage("River Length: " + str(RiverLength))

arcpy.Delete_management(Output_Verticles)
arcpy.Delete_management(Output_Verticles_elevation)
###### END OF CODE TO FIND RIVER SLOPE





# Execute CreateTable (Create the Flood-Discharge Table)
tableName = "FloodDischargeTable.dbf"
arcpy.CreateTable_management(Data_output_Location, tableName)

# Set variables for AddField
out_name = os.path.join(Data_output_Location, tableName)
fieldName1 = "WaterLevel"
fieldName2 = "ChanWidth"
fieldName3 = "CrossArea"
fieldName4 = "WetPerm"
fieldName5 = "HydroRad"
fieldName6 = "RivLength"
fieldName7 = "Q_Min"
fieldName8 = "Q_Med"
fieldName9 = "Q_Maj"

# Execute AddField 9 times for 9 new fields
arcpy.AddField_management(out_name, fieldName1, "DOUBLE", 10, 7, field_is_nullable="NON_NULLABLE")
arcpy.AddField_management(out_name, fieldName2, "DOUBLE", 10, 7, field_is_nullable="NON_NULLABLE")
arcpy.AddField_management(out_name, fieldName3, "DOUBLE", 10, 7, field_is_nullable="NON_NULLABLE")
arcpy.AddField_management(out_name, fieldName4, "DOUBLE", 10, 7, field_is_nullable="NON_NULLABLE")
arcpy.AddField_management(out_name, fieldName5, "DOUBLE", 10, 7, field_is_nullable="NON_NULLABLE")
arcpy.AddField_management(out_name, fieldName6, "DOUBLE", 10, 7, field_is_nullable="NON_NULLABLE")
arcpy.AddField_management(out_name, fieldName7, "DOUBLE", 10, 7, field_is_nullable="NON_NULLABLE")
arcpy.AddField_management(out_name, fieldName8, "DOUBLE", 10, 7, field_is_nullable="NON_NULLABLE")
arcpy.AddField_management(out_name, fieldName9, "DOUBLE", 10, 7, field_is_nullable="NON_NULLABLE")

# Execute DeleteField to remove default field created when CreateTable was executed.
arcpy.DeleteField_management(out_name, ["Field1"])

# Base Discharge Offset variables (set as 0)
Base_Dis = 0


while Water_Level <= MAXwaterLevel:
    
    Water_Level = float(Water_Level + 0.1)

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
    
    # Create Base Offset (original Water_rough was 0.045)
    if Water_Level == 0.1:
        Base_Dis = (1.0 / Water_Rough_flt) * (HydroRadius_flt ** 0.666) * Base_River_Volume * (Slope_flt ** 0.5)

    ###Now to calculate all three of the Discharges
    #Flood Discharge equation (Rough_Q1)
    Q_min = (1.0 / Rough_Q1) * (HydroRadius_flt ** 0.666) * cross_section_area_flt * (Slope_flt ** 0.5) + Base_Dis

    #Flood Discharge equation (Rough_Q2)
    Q_med = (1.0 / Rough_Q2) * (HydroRadius_flt ** 0.666) * cross_section_area_flt * (Slope_flt ** 0.5) + Base_Dis

    #Flood Discharge equation (Rough_Q3)
    Q_maj = (1.0 / Rough_Q3) * (HydroRadius_flt ** 0.666) * cross_section_area_flt * (Slope_flt ** 0.5) + Base_Dis

    #Now to convert all the remaining values to float so that they can be added to the FloodDischarge table
    Q_min_flt = float(Q_min)
    Q_med_flt = float(Q_med)
    Q_maj_flt = float(Q_maj)
    Rivlength_flt = float(RiverLength) #This is here because I didn't know where else to put it



     ###Code to add values to FloodDischarge table


    # Set variables for InsertCursor
    fieldNames = [fieldName1, fieldName2, fieldName3, fieldName4, fieldName5, fieldName6, fieldName7, fieldName8, fieldName9]


    row_values = [(Water_Level, Channel_width_flt, cross_section_area_flt, cross_section_perm_flt,
                                HydroRadius_flt, Rivlength_flt, Q_min_flt, Q_med_flt, Q_maj_flt)]

    # Open an InsertCursor
    cursor_discharge = arcpy.da.InsertCursor(out_name, fieldNames)

    # Insert new rows
    for row in row_values:
        cursor_discharge.insertRow(row)

    # Delete cursor
    del cursor_discharge
    
arcpy.Delete_management(outTable)
arcpy.Delete_management(LandClass_table)