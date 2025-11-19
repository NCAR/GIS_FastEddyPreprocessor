# GIS_FastEddyPreprocessor
These tools will create the inputs for FastEddy using ArcGIS Pro

## Data Requirements
Building Footprints
*	https://gis-fema.hub.arcgis.com/pages/usa-structures
*	https://github.com/microsoft/USBuildingFootprints
*	\\gisData.ucar.edu\data\FacilitiesAndCriticalInfrastructure\buildings

LiDAR
*	https://apps.nationalmap.gov/lidar-explorer/#/
*	Local State, County, or City GIS offices
*	Texas Natural Resources Information Systems - https://tnris.org/stratmap/elevation-lidar/

Elevation
*	\\gisData.ucar.edu\data\ElevationAndDerivedProducts
*	SRTM - \\gisData.ucar.edu\data\ElevationAndDerivedProducts\SRTM\srtm_void_filled\elevation\srtm_n_elev_w.jp2"

Land Cover
*	NLCD - \\gisData.ucar.edu\data\BiologyAndEcology\Landcover_NLCD
*	https://www.mrlc.gov/data

## FastEddy Toolbox - “\customToolbox\FastEddy.pyt"

### Running Tool

1.	Create a new ArcGIS Pro project
2.	Create the following directories in Catalog pane in ArcGIS pro
  * Final
3.	Copy the customTool, lidarTools, all Lidar data into your ArcGIS Pro project
4.	Create the following empty geodatabases:
	* finalGrids
	* scratch
5.	Install Toolbox in ArcGIS Pro
6.	Open ArcGIS Pro and create a new project or open an existing Project
	* In the Catalog Pane, right click on Toolboxes and select Add Toolbox
	* Navigate to the customToolbox directory and select FastEddy.pyt.  Click Ok.

 
Step 1 : Create Building Surface from LiDAR
This tool requires a directory with a bunch of .laz files.  It will unzip the .laz files to .las files and then merge them together into a .lasd file.  Then the tool will create a Digital Terrain Model (DTM), Digital Surface Model (DSM), and a normalized Digital Surface Model (nDSM).


LiDAR EXE Directory from LiDAR Tools 
This is the directory where laszip.exe is located.  This can be downloaded from https://rapidlasso.com/lastools/  This .exe will convert the downloaded .laz files to .las files.	 
Directory where LAS data is stored
This is the directory where you have downloaded and stored the LiDAR .laz files which you have already downloaded.	
Output Geodatabase to store surface data
You must have a geodatabase already created to store your files.  This can be a scratch geodatabase but is where all outputs will be stored.
	
Output LiDAR (lasd) file 
Name of output lasd file and direcotry.  The lasd will be the merged LiDAR files from the input LiDAR directory.	
Input domain name
You specify what the domain name will be for your project.  This will be applied to all outputs. Example Downtown or FortWorth.	

Step 2: Add Height to Buildings
This tool will take the nDSM you created in Step 1 and perform a zonal statistics with the buildings dataset.  The output will be a gridded dataset with the buildings and their heights

a.	Map a Network Connection to //gisdata.ucar.edu/Data
b.	Bring in the Oklahoma Buildings – FacilitiesAndCriticalInfrastructure/buildings/Oklahoma
c.	Select  buildings that fall within the domain and same locally – so you are just working with buildings in Domain – small and faster

Input Surface Data
This is the nDSM that you created in Step 1	 
Buildings with Height output overlaid on top of the LiDAR lasd file.
 
Input Building Height Data
Building footprints that you have already downloaded from FEMA or Microsoft or a local GIS office.	
Input Domain name
This should be the same name you used in Step 1 for consistency.  IE Downtown	
Final Output Directory
This should be a directory that will store the .csv and the domain corner excel file.	
Output Projection File
You may want the projection to match the WRF projection.  You can create a shapefile that is projected to this WRF projection.  Everything will get projected to this projection.	
Scratch Directory
A geodatabase that will be used to store scratch data.	
Domain File
A vector dataset that has the rectangular domain.  Everything will get clipped to this extent.
	
Geodatabase to store final Gridded data
A new geodatabase to store all the clipped and projected data before they get converted to .csv.

 
Step 3: Project, Clip and Convert
This step will take buildings, land cover, land use, and elevation and project it to the WRF projection, clip it to the domain, and convert it to .csv.  Everything will be stored in the final directory you specify

Clip elevation and NLCD to a smaller domain before proceeding to this step.
1.	Zoom into your domain
2.	Bring elevation and NLCD into map
3.	Right click on the layer > Save AS
Input Building Data
The building height raster data that was created in Step 2	 
Input Land Cover Data
Do not use the CONUS level NLCD.  First clip the land cover to a smaller domain and use that as the input.  It will get resampled to 1m, projected, and clipped to the domain.	
Input Elevation Data
Do not use the CONUS level elevation.  First clip the land cover to a smaller domain and use that as the input. 	
Input Land Use Data
Do not use the CONUS level NLCD.  First clip the land cover to a smaller domain and use that as the input. 	
Input domain file
The domain as a vector dataset - used in Step 2	
Input projection file
File that is projected to the WRF projection - used in Step 2	
Final Output Directory
The directory where you want all your .csv files stored. - used in Step 2	
Scratch Output Directory
Scratch geodatabase to store intermittent data	
Geodatabase to store final Gridded Data
The geodatabase to store the projected and clipped gridded data that will eventually be converted to .csv.	








