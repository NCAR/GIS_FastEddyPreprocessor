# FastEddy GIS Prep

## Overview
The purpose of this program is to prepare GIS fields in NetCDF format for use in the FastEddy coupler.

## Instructions
Running the program requires specifying values in the couplerParams.toml parameter file.
1. Domain description. Includes domain name, central lat/lon, width/height, and cell size.
2. Data paths. A base path where all data files reside and the output and log files will be written, as well as names for the elevation, NLCD, and building height files.

Once the parameter file has been filled out, run the program from src via the following command:
{/path/to/python} -m gisPrep {/path/to/parameter/file}

## Notes
Requires Python 3.10+ to read .toml files natively.