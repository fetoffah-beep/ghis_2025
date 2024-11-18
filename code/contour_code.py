# # -*- coding: utf-8 -*-
# """
# Created on Thu Sep 19 14:13:59 2024

# @author: 39351
# """

# # -*- coding: utf-8 -*-
# """
# Created on Thu Sep 19 14:13:59 2024

# @author: 39351
# """

import numpy as np
import matplotlib.pyplot as plt
from osgeo import osr, gdal
from shapely.geometry import LineString
import fiona
from fiona.crs import from_epsg
from pyproj import CRS  # Import pyproj's CRS class

from scipy.interpolate import splprep, splev




# Open the raster
dataset = gdal.Open(r"../DEM/GEE-20240906T150421Z-001/GEE/sample_dem.tif")

# # As described in the Raster Data Model, a GDALDataset contains a list of raster 
# bands, all pertaining to the same area, and having the same resolution. It also 
# has metadata, a coordinate system, a georeferencing transform, size of raster and 
# various other information.


# print("Driver: {}/{}".format(dataset.GetDriver().ShortName,
#                             dataset.GetDriver().LongName))
# print("Size is {} x {} x {}".format(dataset.RasterXSize,
#                                     dataset.RasterYSize,
#                                     dataset.RasterCount))

geotransform = dataset.GetGeoTransform()
# if geotransform:
#     print("Origin = ({}, {})".format(geotransform[0], geotransform[3]))
#     print("Pixel Size = ({}, {})".format(geotransform[1], geotransform[5]))

# # Get raster dimensions
band = dataset.GetRasterBand(1)
# print("Band Type={}".format(gdal.GetDataTypeName(band.DataType)))

srs = osr.SpatialReference()
srs.SetWellKnownGeogCS("WGS84")

# # Example of setting up a dataset (assuming dst_ds is already created)
dataset.SetProjection(srs.ExportToWkt())
# print("Projection is {}".format(dataset.GetProjection()))

# Get geotransformation and projection from the raster
geotransform = dataset.GetGeoTransform()
projection = dataset.GetProjection()


x0, dx, _, y0, _, dy = geotransform
        

band_min = band.GetMinimum()
band_max = band.GetMaximum()
if not band_min or not band_max:
    (band_min, band_max) = band.ComputeRasterMinMax(True)
# print("Min={}, Max={}".format(band_min, band_max))

contour_levels = np.linspace(band_min, band_max, num=50)  # Adjust contour interval

# # print(contour_levels)
        


xsize = band.XSize
ysize = band.YSize

# Define block size to read (e.g., 1024x1024 blocks)
block_size = 10000

plt.figure()

# Schema for the shapefile (line geometries with elevation as property)
schema = {
    'geometry': 'LineString',
    'properties': {'elevation': 'float'}
}

# Define shapefile CRS
crs = CRS.from_epsg(4326).to_wkt()  # WGS84 EPSG code

# Open a new shapefile for writing
with fiona.open('spline_contours.shp', 'w', driver='ESRI Shapefile', crs=crs, schema=schema) as shp:
 

    for i in range(0, xsize, block_size):
        for j in range(0, ysize, block_size):
            # Calculate the block size for the edges
            x_block_size = min(block_size, xsize - i)
            y_block_size = min(block_size, ysize - j)            
            
            # Read the block
            band_array = band.ReadAsArray(i, j, x_block_size, y_block_size)
            
            # Mask invalid data
            band_array = np.ma.masked_invalid(band_array)
            
                    # Skip completely masked blocks
            if band_array.mask.all():
                continue
    
            # Skip blocks with uniform data
            if band_array.min() == band_array.max():
                continue
            
            # # Calculate local min and max for the current block
            # block_min = band_array.min()
            # block_max = band_array.max()
            
            # # Define contour levels based on the block's min and max
            # contour_levels = np.linspace(block_min, block_max, num=50)


            
            # Step 2: Generate X and Y coordinates from the geotransform
            nrows, ncols = band_array.shape
            block_x0 = x0 + dx * i
            block_y0 = y0 + dy * j
            x = np.linspace(block_x0, block_x0 + dx * (x_block_size), x_block_size)
            y = np.linspace(block_y0, block_y0 + dy * (y_block_size), y_block_size)
            X, Y = np.meshgrid(x, y)
    
            # Generate contours
            '''
                to do:  Define the min and max elevations in the global dataset and filter for the ones
                    that fall in the current block           

            '''            
            contour_levels = np.linspace(band_array.min(), band_array.max(), num=50)
            contours=plt.contour(X, Y, band_array, levels=contour_levels)
        
        
        
        
            if not contours.collections:
                print("No contours generated; using default levels.")
                

                
            # Spline-interpolate each contour
            for collection in contours.collections:
                for path in collection.get_paths():
                    # Get the vertices of the contour path
                    vertices = path.vertices
                    # Skip paths with fewer than 4 points
                    if len(vertices) < 4:
                        continue
                    x_vertices, y_vertices = vertices[:, 0], vertices[:, 1]
                    
                    # Apply spline interpolation (tck: spline parameters, u: parametric variable)
                    tck, u = splprep([x_vertices, y_vertices], s=0.0)  # s=0 gives an exact fit
                    u_fine = np.linspace(0, 1, len(u) * 10)  # Create more points for smooth curve
                    x_smooth, y_smooth = splev(u_fine, tck)
            

                    # Create a LineString from the smoothed points
                    smooth_line = LineString(np.column_stack([x_smooth, y_smooth]))

                    # Write the smoothed line to the shapefile with elevation data
                    shp.write({
                        'geometry': smooth_line.__geo_interface__,
                        'properties': {'elevation': float(contours.levels[contours.collections.index(collection)])}
                    })
                    
                    
                                # Plot or store the smooth spline
                    plt.plot(x_smooth, y_smooth)
                    
                    
    
    # Display the combined figure with all blocks
plt.title('Spline-Interpolated Contours for All Blocks')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()

print("Contour shapefile created successfully!")
