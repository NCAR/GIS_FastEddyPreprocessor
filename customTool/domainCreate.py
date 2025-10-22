'''
domainCreate.py

Python code companion to the toolbox tool, domainCreate.pyt.
'''


__date__ = '2024/01/29'
__institution__ = 'NCAR'

import csv
import os

import arcpy


def printMessage(message):
  '''Print message in Geoprocessing progress table.
  
  Parameters
  -----
  messages: message string to be printed

  Returns
  -----
  '''

  arcpy.AddMessage(message)

  return


def createSpatialReference(center_lat, center_lon, domain_name):
  '''Create the FastEddy domain ArcGIS Spatial Reference object.
  The spatial reference will be a Lambert Conformal Conic projection
  using the WGS 1984 (2011) spheroid centered at the specified coordinate.
  
  Parameters
  -----
  center_lat: domain center latitude (degrees N)
  center_lon: domain center longitude (degrees E)
  domain_name: domain name

  Returns
  -----
  ArcGIS SpatialReference
    spatial_ref: FastEddy domain coordinate system as an ArcGIS spatial ref
  '''

# parse spatial reference name
  sr_name = f'FastEddy_{domain_name}_LCC'

# set default parameters
  geo_gcs = 'GCS_WGS_1984'
  datum = 'D_WGS_1984'
  spheroid = 'WGS_1984'
  semimajor = 6378137
  inv_flat = 298.257223563
  proj = 'Lambert_Conformal_Conic'

# parse proj template string
  proj_template = (
    f'PROJCS[{sr_name},'
    f'GEOGCS[{geo_gcs},'
    f'DATUM[{datum},' 
    f'SPHEROID[{spheroid}, {semimajor}, {inv_flat}]],'
    f'PRIMEM["Greenwich", 0.0],'
    f'UNIT["Degree", 0.0174532925199433]],'
    f'PROJECTION[{proj}],'
    f'PARAMETER["False_Easting", 0.0],'
    f'PARAMETER["False_Northing", 0.0],'
    f'PARAMETER["Central_Meridian", {center_lon}],'
    f'PARAMETER["Standard_Parallel_1", {center_lat}],'
    f'PARAMETER["Standard_Parallel_2", {center_lat}],'
    f'PARAMETER["Scale_Factor", 1],'
    f'PARAMETER["Latitude_Of_Origin", {center_lat}],'
    f'UNIT["Meter", 1.0]]'
  )

  spatial_ref = arcpy.SpatialReference(text=proj_template)

  return (spatial_ref)


def domainCenter(domain_name, out_gdb, spatial_ref):
  '''Creates a domain center point as a new point feature class.
  
  Parameters
  -----
  domain_name: domain name
  out_gdb: path to output geodatabase
  spatial_ref: domain spatial reference, e.g. output by createSpatialReference
  
  Returns
  -----
  ArcGIS Point Feature Class
    domain_center_fc: domain center point as an ArcGIS feature class
  '''

# parse feature class name
  domain_center_name = f'{domain_name}_center'

# create empty point feature class
  domain_center_fc = arcpy.management.CreateFeatureclass(
    out_gdb, domain_center_name, 
    'POINT', '', 'DISABLED', 'DISABLED', 
    spatial_ref
  )

# insert domain center point into feature class
  with arcpy.da.InsertCursor(domain_center_fc, ['SHAPE@']) as cursor:
    domain_center_point = arcpy.PointGeometry(
      arcpy.Point(), spatial_ref
    )
    cursor.insertRow((domain_center_point))
  
  return (domain_center_fc)


def domainPolygon(
    domain_center_fc, domain_width, domain_height, domain_name, out_gdb
  ):
  '''Calculates the domain bounds in the form of a polygon feature class.
  Step 1: create buffer.
  Step 2: wrap buffer with polygon.
  Works for square domains only; rectangles require different logic.

  Parameters
  -----
  domain_center_fc: point feature class containing the projected domain center
  domain_width: width of domain (km)
  domain_height: height of domain (km)
  domain_name: domain name
  out_gdb: path to output geodatabase

  Returns
  -----
  ArcGIS Polygon Feature Class
    domain_polygon: domain as a polygon feature class
  '''

# parse paths to domain buffer and feature envelope
  domain_buffer_path = os.path.join(out_gdb, f'{domain_name}_FE_buffer')
  domain_path = os.path.join(out_gdb, f'{domain_name}_FE_domain')

# create buffer
  domain_buffer = arcpy.analysis.Buffer(
    domain_center_fc, domain_buffer_path, domain_width
  )

# create feature envelope around domain buffer
  domain_polygon = arcpy.management.FeatureEnvelopeToPolygon(
    domain_buffer, domain_path
  )

  return (domain_polygon)


def domainVertices(domain_polygon, domain_name, out_gdb, spatial_ref, save):
  '''Calculates domain vertice lat/lon coordinates.
  Optionally, save select output to csv.
  
  Parameters
  -----
  domain_polygon: polygon feature class, e.g. output of domainPolygon
  domain_name: domain name
  out_gdb: path to output geodatabase
  spatial_ref: domain spatial reference, e.g. output by createSpatialReference
  save: option to save select output to csv
  
  Returns
  -----
  ArcGIS Point Feature Class
  '''

# parse path to vertices feature class
  vertices_path = os.path.join(out_gdb, f'{domain_name}_vertices')

# calculate domain vertices
  domain_vertices = arcpy.management.FeatureVerticesToPoints(
    domain_polygon, vertices_path
  )

# add xy coordinates
  arcpy.management.AddXY(domain_vertices)

# add empty lat/lon fields
  arcpy.management.AddFields(
    domain_vertices, [
      ['Lat', 'DOUBLE', 'Latitude'],
      ['Lon', 'DOUBLE', 'Longitude']
    ]
  )

# populate lat/lon fields
  arcpy.management.CalculateGeometryAttributes(
    in_features=domain_vertices,
    geometry_property=[['lat', 'POINT_Y'], ['lon', 'POINT_X']],
    coordinate_system=spatial_ref,
    coordinate_format='DD'
  )

# add and populate Name field
  arcpy.management.AddField(
    in_table=domain_vertices,
    field_name='Name',
    field_type='TEXT'
  )
  with arcpy.da.UpdateCursor(domain_vertices, 'Name') as cursor:
    for row in cursor:
      row[0] = domain_name
      cursor.updateRow(row)

# if save is True, save output to csv
  if bool(save) == True:

# save select output to csv
    keep_fields = ['POINT_X', 'POINT_Y', 'Lat', 'Lon', 'Name']
    data = [keep_fields]
    with arcpy.da.SearchCursor(domain_vertices, keep_fields) as cursor:
      for row in cursor:
        data.append(list(row))

    with open(f'DomainBoundaries_{domain_name}.csv', 'w', newline='') as f:
      writer = csv.writer(f)
      writer.writerows(data)

  return (domain_vertices)