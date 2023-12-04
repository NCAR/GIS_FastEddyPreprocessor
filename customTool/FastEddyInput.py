
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# ** Copyright UCAR (c) 2018
# ** University Corporation for Atmospheric Research(UCAR)
# ** National Center for Atmospheric Research(NCAR)
# ** Research Applications Laboratory(RAL)
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
# ** 2016/4/27 10:00:00
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

# --- Import Modules --- #
import sys
import os
import csv                                                                      # For reading input forecast points in CSV format
import time
import math
import zipfile
from zipfile import ZipFile, ZipInfo
import shutil
from collections import defaultdict                                             # Added 09/03/2015 Needed for topological sorting algorthm
from itertools import takewhile, count                                          # Added 09/03/2015 Needed for topological sorting algorthm                                                                  # Added 07/10/2015 for ExmineOutputs
import netCDF4
import numpy
from arcpy.sa import *
import subprocess
import shutil
import glob
import pandas as pd

# --- End Import Modules --- #




#variables used in the creation of surface data from lidar data
cell_size = 1
only_ground_plus_class_code = False
class_code = 6
minimum_height = -5
maximum_height = 3000
processing_extent = True
classify_noise = False
log_directory = True
verbose = True

fullExtent = ""
snap = ""

def convertlas(arcpy,laszip,inputLidarDir,lasdName):
        
        printMessages(arcpy, ['Step 1: Convert laz to las...'])
        cmd = "%s %s\*.laz"%(laszip,inputLidarDir)
        ExecuteProcess(arcpy,cmd)
        
        #lasfiles = [f for f in os.listdir(inputLidarDir) if f.endswith('.las')]
        lasfiles = sorted(glob.glob(f'{inputLidarDir}\\*.las'))
        arcpy.env.workspace = inputLidarDir
        printMessages(arcpy, ['Step 2: Create lasd'])
        arcpy.management.CreateLasDataset(lasfiles, lasdName)
        out_las_dataset = os.path.join(inputLidarDir,lasdName)
        
        
        return out_las_dataset


def get_name_from_feature_class(arcpy,feature_class):
    desc_fc = arcpy.Describe(feature_class)
    return desc_fc.name

def create_msg_body(msg_prefix, start_time, end_time):

    # Creates the message returned after each run of a function (successful or unsuccessful)
    diff = end_time - start_time

    if diff > 0:
        if msg_prefix == "":
            msg_prefix = "Elapsed time: "
        else:
            msg_prefix = msg_prefix + "  Elapsed time: "

        elapsed_time_mins = int(math.floor(diff/60))
        minutes_txt = " minutes "
        if elapsed_time_mins == 1:
            minutes_txt = " minute "
        if elapsed_time_mins > 0:
            elapsed_time_secs = int(round(diff - (60 * elapsed_time_mins)))
            seconds_txt = " seconds."
            if elapsed_time_secs == 1:
                seconds_txt = " second."
            elapsed_time_formatted = str(elapsed_time_mins) + minutes_txt + str(elapsed_time_secs) + seconds_txt
        else:
            elapsed_time_secs = round(diff - (60 * elapsed_time_mins), 2)
            seconds_txt = " seconds."
            if elapsed_time_secs == 1:
                seconds_txt = " second."
            elapsed_time_formatted = str(elapsed_time_secs) + seconds_txt

        msg_body = msg_prefix + elapsed_time_formatted

    else:
        msg_body = msg_prefix

    return msg_body


def extract(arcpy,lc_lasd, lc_ws, lc_cell_size, lc_ground_classcode,  lc_class_code, lc_output_elevation, lc_minimum_height,
            lc_maximum_height, lc_processing_extent, lc_noise, lc_log_dir, lc_debug, lc_memory_switch):

    dem = None
    dsm = None
    ndsm = None

           

    # create dem
    desc = arcpy.Describe(lc_lasd)
    l_unit = desc.spatialReference.linearUnitName
   
    #        if desc.spatialReference.linearUnitName in ['Foot_US', 'Foot']:
    if 'feet' in l_unit.lower() or 'foot' in l_unit.lower():
        unit = 'Feet'
    else:
        unit = 'Meters'

    if lc_class_code == 15:
        lc_cell_size = lc_cell_size*2

    # Classify overlap points
    # ptSpacing = desc.pointSpacing * 2.25
    # sampling = '{0} {1}'.format(ptSpacing, unit)
    # arcpy.ClassifyLasOverlap_3d(lc_lasd, sampling)

    # get lidar class code - TEMPORARY until Pro 2.3
    msg_body = create_msg_body("Looking for class codes: ", 0, 0)
    print(msg_body)

    las_desc = arcpy.Describe(lc_lasd)
    class_codes = las_desc.classCodes
    class_code_list = [int(code) for code in class_codes.split(';')]

    ground_code = 2

    if lc_ground_classcode and lc_class_code in class_code_list:
        class_code_list = list()
        class_code_list.append(int(ground_code))
        class_code_list.append(int(lc_class_code))

    # Generate DEM
    if ground_code in class_code_list:
        dem = lc_output_elevation + "_dtm"

        if arcpy.Exists(dem):
            arcpy.Delete_management(dem)

        msg_body = create_msg_body("Creating Ground Elevation using the following class codes: " +
                                       str(ground_code), 0, 0)
        print(msg_body)

        ground_ld_layer = arcpy.CreateUniqueName('ground_ld_lyr')

        # Filter for ground points
        arcpy.management.MakeLasDatasetLayer(lc_lasd, ground_ld_layer, class_code=str(ground_code))

        arcpy.conversion.LasDatasetToRaster(ground_ld_layer, dem, 'ELEVATION',
                                                'BINNING MAXIMUM LINEAR',
                                                sampling_type='CELLSIZE',
                                                sampling_value=lc_cell_size)

        lc_max_neighbors = "#"
        lc_step_width = "#"
        lc_step_height = "#"

        if lc_noise:
            # Classify noise points
            msg_body = create_msg_body("Classifying points that are " + lc_minimum_height + " below ground and " +
                                           lc_maximum_height + " above ground as noise.", 0, 0)
            print(msg_body)

            arcpy.ClassifyLasNoise_3d(lc_lasd, method='RELATIVE_HEIGHT', edit_las='CLASSIFY',
                                          withheld='WITHHELD', ground=dem,
                                          low_z=lc_minimum_height, high_z=lc_maximum_height,
                                          max_neighbors=lc_max_neighbors, step_width=lc_step_width,
                                          step_height=lc_step_height,
                                          extent=lc_processing_extent)
        else:
            # Classify noise points
            printMessages(arcpy,["Noise will not be classified.", 0, 0])
            

        # check if we need to create dsm and ndsm based on lc_class_code != -1

        if lc_class_code != -1:
            # create dsm
            dsm = lc_output_elevation + "_dsm"

            if arcpy.Exists(dsm):
                arcpy.Delete_management(dsm)

            printMessages(arcpy,["Creating Surface Elevation using the following class codes: " +
                                           str(class_code_list), 0, 0])
            

            dsm_ld_layer = arcpy.CreateUniqueName('dsm_ld_lyr')

            return_usage = arcpy.Usage("MakeLasDatasetLayer_management").split(', ')[3].strip('{}').split(' | ')
            # last return = first entry
            last_return = return_usage[0]

            if lc_class_code == 15:
                arcpy.management.MakeLasDatasetLayer(lc_lasd, dsm_ld_layer, class_code=class_code_list)
            else:
                arcpy.management.MakeLasDatasetLayer(lc_lasd, dsm_ld_layer, class_code=class_code_list,
                                                         return_values=[last_return])

            arcpy.conversion.LasDatasetToRaster(in_las_dataset=dsm_ld_layer,
                                                    out_raster=dsm,
                                                    value_field='ELEVATION',
                                                    interpolation_type='BINNING MAXIMUM LINEAR',
                                                    sampling_type='CELLSIZE',
                                                    sampling_value=lc_cell_size)

            # create ndsm
            printMessages(arcpy,["Creating normalized Surface Elevation using " +
                                           get_name_from_feature_class(arcpy,dsm) + " and " +
                                           get_name_from_feature_class(arcpy,dem), 0,0])
            

            ndsm = lc_output_elevation + "_ndsm"

            if arcpy.Exists(ndsm):
                arcpy.Delete_management(ndsm)

            arcpy.Minus_3d(dsm, dem, ndsm)
        else:
            printMessages(arcpy,["Couldn't detect ground class code in las dataset. "
                                       "Use the -Classify LAS Ground- tool to classify ground. "
                                       "Exiting...", 0, 0])
            

    return dem, dsm, ndsm

def buildingHeight():
        #Buildings
        structure_name_proj = "structure_project_%s"%domainName
        structure_name = "structure_%s"%domainName
        out_structure_ascii_name = "building_%s.asc"%domainName
        projbuilding = os.path.join(scratchDir,structure_name_proj)
        out_building = os.path.join(WRFDataDir,structure_name)
        out_building_ascii =os.path.join(finalDir,out_structure_ascii_name)

        arcpy.management.Project(buildings, projbuilding, spatial_ref)
        #arcpy.analysis.Clip(projbuilding, out_domain, out_building)

        print("created %s"%projbuilding)
        #Zonal with buildings
        arcpy.env.cellAlignment = "ALIGN_WITH_PROCESSING_EXTENT"
        with arcpy.EnvManager(extent=fullExtent, snapRaster=out_height):
            domain_building = arcpy.sa.ZonalStatistics(in_zone_data=projbuilding, zone_field="OBJECTID", in_value_raster=projheight, statistics_type="MAXIMUM", ignore_nodata="DATA", process_as_multidimensional="CURRENT_SLICE", percentile_value=90, percentile_interpolation_type="AUTO_DETECT")
            domain_building.save(out_building)
            
        arcpy.conversion.RasterToASCII(out_building, out_building_ascii)
        print("finished with buildings %s"%out_building_ascii)

def projection(arcpy,inData,scratchDir,projection,domainName):
        #Make sure the domain and all data are projected to WRF Lambert       

        spatial_ref = arcpy.Describe(projection).spatialReference
        buildingsBase = os.path.basename(inData)
        name = buildingsBase.split(".")[0]
        outName = name + "_"+domainName 
        ###make sure building 
        out_data = os.path.join(scratchDir,outName)
        #project the domain to WRF Lambert
        printMessages(arcpy,["projectiong vector data " + inData
                                            + " to become " +
                                           out_data, 0,0])
        
        arcpy.management.Project(inData, out_data, spatial_ref)

        return out_data

def clipRaster(arcpy,finalDir,domainRaster,inData,scratchDir,outDir,domain,domainName,dataName):
        #Make sure the domain and all data are projected to WRF Lambert       
        dataBase = os.path.basename(inData)
        name = dataBase.split(".")[0]
        outName = dataName + "_" + domainName
        outName2 = outName + "_resample"
        ###make sure building
        out_raster = os.path.join(scratchDir,outName2)
        out_data = os.path.join(outDir,outName)
        out_ascii = os.path.join(finalDir,outName+".asc")
        
        #project the domain to WRF Lambert
        printMessages(arcpy,["resample raster data " + inData
                                            + " and " +
                                           out_raster, 0,0])
        ##resample
        arcpy.env.cellAlignment = "ALIGN_WITH_PROCESSING_EXTENT"
        with arcpy.EnvManager(snapRaster=domainRaster):
                arcpy.management.Resample(inData, out_raster, 1)
                #resampled_raster = arcpy.sa.Resample(inData, "NearestNeighbor", "NONE", 1)
                #resampled_raster.save =(out_raster)
                
        printMessages(arcpy,["clipping raster data " + out_raster
                                            + " and " +
                                           out_data, 0,0])

        with arcpy.EnvManager(extent=domainRaster, snapRaster=domainRaster):
                #arcpy.management.Clip(out_raster, domain, out_data)
                mask = arcpy.sa.ExtractByMask(in_raster=out_raster, in_mask_data=domainRaster, extraction_area="INSIDE")
                mask.save(out_data)
                
        printMessages(arcpy,["converting to ascii " + out_data + " and " + out_ascii, 0,0])

        with arcpy.EnvManager(extent=domainRaster, snapRaster=domainRaster):
                arcpy.conversion.RasterToASCII(out_data, out_ascii)
                
        return out_data

def clipBuildings(arcpy,inData,scratchDir,domain):
        #Make sure the domain and all data are projected to WRF Lambert       
        dataBase = os.path.basename(inData)
        name = dataBase.split(".")[0]
        outName = name + "_clip"
        ###make sure building 
        out_data = os.path.join(scratchDir,outName)
        #project the domain to WRF Lambert
        printMessages(arcpy,["clipping data " + inData
                                            + " to become " +
                                           out_data, 0,0])
        ##Add clip
        arcpy.analysis.Clip(inData, domain, out_data)
        return out_data


def projectionRaster(arcpy,inData,scratchDir,projection,domainName,dataName):
        #Make sure the domain and all data are projected to WRF Lambert       

        spatial_ref = arcpy.Describe(projection).spatialReference
        
        outName = dataName + "_" + domainName + "_rasterProject"
        ###make sure building 
        out_data = os.path.join(scratchDir,outName)
        #project the domain to WRF Lambert
        printMessages(arcpy,["projectiong raster data " + inData
                                            + " and " +
                                           out_data, 0,0])
        
        arcpy.env.cellAlignment = "ALIGN_WITH_PROCESSING_EXTENT"
        with arcpy.EnvManager(cellSize=1, extent=fullExtent, snapRaster=snap):
            raster = arcpy.management.ProjectRaster(inData, out_data, spatial_ref)
            
        
        return out_data

def getHeights(arcpy,building,surface,dataOutDir,domainName,domain ):
        #get heights for building
        with arcpy.da.SearchCursor(domain, ["SHAPE@"]) as cursor:
                for row in cursor:
                        fullExt = row[0].extent

        fullExtent = fullExt
        snap = surface
        height_name_proj = "building_"+ domainName   
        outputData = os.path.join(dataOutDir, height_name_proj)
        printMessages(arcpy,["creating " + outputData, 0, 0])
        
        arcpy.env.cellAlignment = "ALIGN_WITH_PROCESSING_EXTENT"
        with arcpy.EnvManager(cellSize=1,extent=fullExtent, snapRaster=snap):
            domain_building = arcpy.sa.ZonalStatistics(in_zone_data=building, zone_field="OBJECTID", in_value_raster=surface, statistics_type="MAXIMUM", ignore_nodata="DATA", process_as_multidimensional="CURRENT_SLICE", percentile_value=90, percentile_interpolation_type="AUTO_DETECT")
            domain_building.save(outputData)
       

def domainCorners(arcpy, domainName, scratchDir, finalDir, projection, out_domain):
        #Make sure the domain and all data are projected to WRF Lambert

# DOMAIN information
        cornerPT = "corners_%s"%domainName
        domain_pt = os.path.join(scratchDir, cornerPT)
        excelname = "DomainBoundaries_%s.xlsx"%domainName
        domain_excel = os.path.join(finalDir, excelname)
        csvname = "DomainBoundaries_%s.csv"%domainName
        domain_csv = os.path.join(finalDir, csvname)
# status message        
        with arcpy.da.SearchCursor(out_domain, ["SHAPE@"]) as cursor:
                for row in cursor:
                    fullExtent = row[0].extent
        printMessages(arcpy,["getting full extent " +str(fullExtent), 0,0])
# create domain vertices
        arcpy.management.FeatureVerticesToPoints(
            in_features=out_domain, 
            out_feature_class=domain_pt, 
            point_location="ALL"
        )
# calculate coordinates for domain vertices
        spatial_ref = arcpy.da.Describe(projection)['spatialReference']
        arcpy.management.AddXY(domain_pt)
        arcpy.management.AddFields(domain_pt, [
            ['Lat', 'DOUBLE', 'Latitude'],
            ['Lon', 'DOUBLE', 'Longitude']
        ])
        arcpy.management.CalculateGeometryAttributes(
            in_features=domain_pt,
            geometry_property=[['lat', 'POINT_Y'], ['lon', 'POINT_X']],
            coordinate_system=spatial_ref,
            coordinate_format='DD'
        )
# add and populate Name field
        arcpy.management.AddField(
            in_table=domain_pt,
            field_name='Name',
            field_type='TEXT'
        )
        with arcpy.da.UpdateCursor(domain_pt, 'Name') as cursor:
            for row in cursor:
                row[0] = domainName
                cursor.updateRow(row)  
# save output to excel
        arcpy.conversion.TableToExcel(
            Input_Table=domain_pt, 
            Output_Excel_File=domain_excel, 
            Use_field_alias_as_column_header="NAME", 
            Use_domain_and_subtype_description="CODE"
        )
# remove unwanted fields from table
        keep_fields = ['POINT_X', 'POINT_Y', 'Lat', 'Lon', 'Name']
        domain_pd = pd.read_excel(domain_excel)
        drop_fields = [x for x in domain_pd.columns if x not in keep_fields]
        domain_pd.drop(drop_fields, axis=1, inplace=True)
        domain_pd.to_csv(domain_csv, index=False)
        
##def domainCorners(arcpy, domainName,scratchDir, finalDir,projection,out_domain):
##        #Make sure the domain and all data are projected to WRF Lambert
##
##        #DOMAIN information
##        cornerPT = "corners_%s"%domainName
##        domain_pt = os.path.join(scratchDir, cornerPT)
##        excelname = "DomainBoundaries_%s.xlsx"%domainName
##        domain_excel= os.path.join(finalDir, excelname)
##        
##        spatial_ref = arcpy.Describe(projection).spatialReference
##        
##        with arcpy.da.SearchCursor(out_domain, ["SHAPE@"]) as cursor:
##                for row in cursor:
##                    fullExtent = row[0].extent
##        printMessages(arcpy,["getting full extent " +str(fullExtent), 0,0])
##
##        arcpy.management.FeatureVerticesToPoints(in_features=out_domain, out_feature_class=domain_pt, point_location="ALL")
##
##        # Generate the extent coordinates using Add Geometry Properties tool
##        coord_sys = arcpy.SpatialReference(4326) # EPSG:4326 -> GCS_WGS_1984
##        #arcpy.AddGeometryAttributes_management(domain_pt, "POINT_X_Y_Z_M", "", "", coord_sys)
##        arcpy.management.AddFields(domain_pt, [['name','TEXT','name']])
##        
##        # Process: Calculate Geometry Attributes (Calculate Geometry Attributes) (management)
##        corners_test = arcpy.management.CalculateGeometryAttributes(in_features=domain_pt, geometry_property=[["lat", "POINT_Y"], ["lon", "POINT_X"]], length_unit="", area_unit="", coordinate_system="GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]]", coordinate_format="DD")[0]
##        with arcpy.da.UpdateCursor(domain_pt, ["name"]) as cursor:
##            for row in cursor:
##                row[0]  = domainName
##        #        row[2] = row[1] #cal lat to be Point_Y
##        #        row[4] = row[3]#cal lon to be Point_X
##                cursor.updateRow(row)
##                
##        # Process: Add XY Coordinates (Add XY Coordinates) (management)
##        arcpy.management.AddXY(in_features=domain_pt)[0]
##        arcpy.management.FeatureVerticesToPoints(in_features=out_domain, out_feature_class=domain_pt, point_location="ALL")
##
##
####        # Generate the extent coordinates using Add Geometry Properties tool
####        coord_sys = arcpy.SpatialReference(4326) # EPSG:4326 -> GCS_WGS_1984
####        arcpy.AddGeometryAttributes_management(domain_pt, "POINT_X_Y_Z_M", "", "", coord_sys)
####        arcpy.management.AddFields(domain_pt, [['lat', 'DOUBLE', 'lat'], ['lon', 'DOUBLE', 'lon']])
####        with arcpy.da.UpdateCursor(domain_pt, ["name","POINT_Y","lat","POINT_X","lon"]) as cursor:
####            for row in cursor:
####                row[0]  = domainName
####                row[2] = row[1] #cal lat to be Point_Y
####                row[4] = row[3]#cal lon to be Point_X
####                print(row[1])
####                cursor.updateRow(row)
####                
####        # Process: Add XY Coordinates (Add XY Coordinates) (management)
####        arcpy.management.AddXY(in_features=domain_pt)[0]
##        # Process: Table To Excel (Table To Excel) (conversion)
##        arcpy.conversion.TableToExcel(Input_Table=domain_pt, Output_Excel_File=domain_excel, Use_field_alias_as_column_header="NAME", Use_domain_and_subtype_description="CODE")

        
def convert2excel(inData,outTable):
        printMessages(arcpy,["converting data to excel file name", 0,0])
        arcpy.conversion.TableToExcel(Input_Table=domain_pt, Output_Excel_File=domain_excel, Use_field_alias_as_column_header="NAME", Use_domain_and_subtype_description="CODE")
        

def printMessages(arcpy, messages):
    '''provide a list of messages and they will be printed and returned as tool messages.'''
    #for i in messages:
    #    print(i)
    #    arcpy.AddMessage(i)
    arcpy.AddMessage(messages)
        
def ExecuteProcess(arcpy,string):
    p = subprocess.Popen(string, shell=True)
    ph_ret = p.wait()
    
    

