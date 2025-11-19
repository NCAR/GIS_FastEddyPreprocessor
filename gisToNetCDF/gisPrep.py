'''
gisPrep.py

Prepares a GIS NetCDF file for ingest by the FastEddy coupler.
Uses data output by ArcGIS Pro.
'''

import logging
from pathlib import Path
import sys
import tomllib

import numpy as np
from pyproj import CRS, Transformer
import xarray as xr


class GISPrep:

  def __init__(self, pf):
    '''Loads in and parses parameters and creates projection.'''

    # load parameters
    with open(pf, 'rb') as pf:
      self.params = tomllib.load(pf)

    # convert base_path to pathlib.Path object
    base_path = Path(self.params['base_path'])
    
    # initialize log file
    log_file_name = f'{self.params['domain_name']}_log.txt'
    log_file = base_path.joinpath(log_file_name)
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
      )
    
    # build input and output file paths
    self.elev_fp = base_path.joinpath(self.params['elev_file'])
    logging.info(f'Elevation data file: {self.elev_fp}')
    self.nlcd_fp = base_path.joinpath(self.params['nlcd_file'])
    logging.info(f'NLCD data file: {self.nlcd_fp}')
    self.buildings_fp = base_path.joinpath(self.params['buildings_file'])
    logging.info(f'Building height data file: {self.buildings_fp}')
    self.nc_fp = base_path.joinpath(f'{self.params['domain_name']}_gis.nc')
    logging.info(f'NetCDF output file: {self.nc_fp}')

    # parse proj-string
    self.proj_string = (
      f'+proj=lcc +lon_0={self.params['lon_0']} +lat_0={self.params['lat_0']} '
      f'+lat_1={self.params['lat_0']} +lat_2={self.params['lat_0']}'
    )
    logging.info(f'proj-string: {self.proj_string}')

    return


  def genCoords(self):
    '''Generates X, Y, LON, and LAT coordinate arrays.'''

    # create x and y indices
    logging.info('Creating x and y indices...')
    self.xs = np.arange(
      -self.params['domain_width']/2, self.params['domain_width']/2
    )
    self.ys = np.arange(
      -self.params['domain_height']/2, self.params['domain_height']/2
    )
    xy = np.meshgrid(self.xs, self.ys)
    logging.info('Index creation complete.')

    # transform to geodetic coordinates
    logging.info('Transforming to WGS 1984 geodetic coordinates...')
    crs = CRS.from_proj4(self.proj_string)
    t = Transformer.from_crs(crs, 4326)  # 4326: WGS 1984
    self.lats, self.lons = t.transform(xy[0], xy[1])
    logging.info('Coordinate transform complete.')

    return


  def gisNetCDF(self):
    '''Creates GIS NetCDF.'''

    # open input Datasets
    elev = np.loadtxt(self.elev_fp, skiprows=6)
    nlcd = np.loadtxt(self.nlcd_fp, skiprows=6)
    buildings = np.loadtxt(self.buildings_fp, skiprows=6)

    # flip GIS fields along the y-axis (move origin from NW to SW)
    elev = np.flip(elev, axis=0)
    nlcd = np.flip(nlcd, axis=0)
    buildings = np.flip(buildings, axis=0)

    # build Dataset
    ds = xr.Dataset(
      data_vars=dict(
        x=(['x'], self.xs),
        y=(['y'], self.ys),
        elevation=(['y', 'x'], elev),
        lat=(['y', 'x'], self.lats),
        lon=(['y', 'x'], self.lons),
        LandCover=(['y', 'x'], nlcd),
        BuildingHeights=(['y', 'x'], buildings),
        cellsize=self.params['cellsize']
      ),
      attrs=dict(
        description='NetCDF file created from ArcGIS inputs',
      )
    )

    # Save Dataset
    ds.to_netcdf(self.nc_fp)

    logging.info('Dataset creation complete. Good bye.')

    return


if __name__ == '__main__':
  if len(sys.argv) == 2:
    GP = GISPrep(sys.argv[1])
    logging.info(f'Preparing FastEddy GIS for {GP.params['domain_name']}')
  else:
    print(
      'Usage: python -m validate {path_to_parameter_file}'
    )
    exit()
  
  GP.genCoords()
  GP.gisNetCDF()