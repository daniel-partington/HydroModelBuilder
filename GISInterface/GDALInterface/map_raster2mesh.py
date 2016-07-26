# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 08:11:11 2016

@author: part0075

Mesh builder from stratigraphy

Build from bottom to top.

This is only to be called once the rasters have been reprojected into the
correct coordinate system

"""

from osgeo import gdal, gdalconst
import array2Vtk
#import imp
#imp.load_source('build_vtk_from_array', '../../Utilities/array2Vtk.py')

# Step 1. Load raster layers top and bottom

#import matplotlib.pyplot as plt

def map_raster_array_to_mesh(hu_raster_path, hu_raster_files, out_path, vtk_out, min_height, max_height):
    raster_set = {}
    
    for raster in hu_raster_files:
        fname = hu_raster_path + raster # + hu_ext
        #print 'Processing: ', fname    
        ds = gdal.Open(fname, gdalconst.GA_ReadOnly)    
        raster_set[raster] = [ds.GetRasterBand(1).ReadAsArray(), 0.0] 
        #raster_set[raster] = [ds.GetRasterBand(1).ReadAsArray(), ds.GetRasterBand(1).GetNoDataValue()]
        #plt.figure()
        #plt.imshow(raster_set[raster][0])            

        if raster == hu_raster_files[0]:
            ncol = ds.RasterXSize
            nrow = ds.RasterYSize       
            (x0, delr, z0, y0, z1, delc) = ds.GetGeoTransform()
            delc = -delc
        ds = None
    
    # raster_extrap
    import numpy as np
    from scipy import ndimage as nd
    
    # This solution taken from: http://stackoverflow.com/questions/3662361/fill-in-missing-values-with-nearest-neighbour-in-python-numpy-masked-arrays 
    def fill(data, invalid=None):
        """
        Replace the value of invalid 'data' cells (indicated by 'invalid') 
        by the value of the nearest valid data cell
    
        Input:
            data:    numpy array of any dimension
            invalid: a binary array of same shape as 'data'. True cells set where data
                     value should be replaced.
                     If None (default), use: invalid  = np.isnan(data)
    
        Output: 
            Return a filled array. 
        """
    
        if invalid is None: invalid = np.isnan(data)
    
        ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
        return data[tuple(ind)]
    
    #raster_statistics = {}
    
    for raster in reversed(hu_raster_files):
        data = raster_set[raster][0]
        invalid = np.ma.masked_array(data, data==raster_set[raster][1])
        raster_set[raster] += [invalid]
        a = fill(data , invalid=invalid.mask)
        raster_set[raster] += [a]
        #raster_statistics[raster] = [np.min(data[~invalid.mask]), np.max(data[~invalid.mask]), np.mean(data[~invalid.mask])]    
        #Clean_up
        data = None
        invalid = None
        a = None
    
    
    layers = len(hu_raster_files)/2
    
    mesh = np.zeros((layers+1, 
                     len(raster_set[hu_raster_files[0]][0]), 
                     len(raster_set[hu_raster_files[0]][0][0])))
    
    
    zone_matrix = np.zeros((layers, 
                          len(raster_set[hu_raster_files[0]][0]), 
                          len(raster_set[hu_raster_files[0]][0][0])))
    
    thickness = np.zeros((layers, 
                          len(raster_set[hu_raster_files[0]][0]), 
                          len(raster_set[hu_raster_files[0]][0][0])))
                          
    # Start by setting bottom 
                          
    #thickness_statistics ={}
    min_thick_mask = []    
    
    for index, raster in enumerate(reversed(hu_raster_files)):
        # Only work on even index
        if index%2 == 1:
            continue
        
        if index == 0: # This is just for setting the bottom
            # Set base        
            mesh[0] = raster_set[raster][3]
            # Set top of bottom layer     
            mesh[1] = raster_set[hu_raster_files[len(hu_raster_files)-index-2]][3]

            # Check for minimum thickness and adjust if needed
            mask_active = raster_set[raster][2].mask
            mask_thickness = mesh[1]-mesh[0] < min_height
            
            mask_adjust = mask_active & mask_thickness              
            
            min_thick_mask += [mask_thickness]
            
            mesh[1][mask_thickness] = mesh[0][mask_thickness] + min_height
        else:        
            # Check bottom of next layer against top of previous layer        
            
            # For values not already defined in top of previous layer modify elevations
            mask = raster_set[hu_raster_files[len(hu_raster_files)-index+1]][2].mask
            mask2 = raster_set[raster][2].mask
            mask3 = mask & mask2
            #mask3 = ~mask3                

            # Check if when changing value for bottom of current layer that it
            # isn't pushing the top of the layer below its bottom.
            push_below_mask = mesh[index/2] > raster_set[raster][3]
            push_below_mask_in_context = push_below_mask & mask
            #print push_below_mask_in_context
            mesh[index/2][mask] = raster_set[raster][3][mask]

            #mesh[index/2][~mask3] = raster_set[raster][3][~mask3]
            mesh[index/2+1] = raster_set[hu_raster_files[len(hu_raster_files)-index-2]][3] # hu_raster_files[index+1]
            # Fix up for where top is below bottom for layer
            
            mesh[index/2+1][mesh[index/2+1] < mesh[index/2]] = mesh[index/2][mesh[index/2+1] < mesh[index/2]] + min_height

            # Check for minimum thickness and adjust if needed
            mask_thickness = mesh[index/2+1]-mesh[index/2] < min_height
            min_thick_mask += [mask_thickness]
            
            mask_active = raster_set[raster][2].mask

            mask_adjust = mask_active & ~mask_thickness              

            mesh[index/2+1][mask_thickness] = mesh[index/2][mask_thickness] + min_height
        #End if
            
        zone_matrix[index/2][raster_set[raster][2].mask] = -1        
        zone_matrix[index/2][~raster_set[raster][2].mask] = index/2+1        
    #End for
        
    import matplotlib.pyplot as plt
    for index, raster in enumerate(reversed(hu_raster_files)):
        # Only work on even index
        #print index, raster    
        if index%2 == 1:
            continue
    
        thickness[index/2] = mesh[index/2+1]-mesh[index/2]        
        #thickness_statistics[raster] = [thickness[index/2], np.mean(thickness[index/2][~raster_set[raster][2].mask])]
    # End for
        
    grid_info = [ncol, nrow, delc, delr, x0, y0]
    array2Vtk.build_vtk_from_array(grid_info, np.fliplr(mesh), ["z_elev"], [np.fliplr(mesh)], ["zone", "thickness"], [np.fliplr(zone_matrix), np.fliplr(thickness)], out_path, vtk_out)

    mesh = np.flipud(mesh)
    zone_matrix = np.flipud(zone_matrix)

    return mesh, zone_matrix


if __name__ == "__main__":
    bottom_ext = "_2b"
    top_ext = "_1t"
    
    #hu_raster_path = r"C:\Workspace\part0075\MDB modelling\VAF_v2.0_ESRI_GRID\ESRI_GRID\Preprocessed_data\Hydrogeological_Unit_Layers\\"
    #hu_raster_files = ["qa_1t_bb", "qa_2b_bb", "utqa_1t_bb", "utqa_2b_bb", "utaf_1t_bb", "utaf_2b_bb", "lta_1t_bb", "lta_2b_bb", "cps_1t_bb", "cps_2b_bb"]
    
    hu_raster_path = r"C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\model_files\\" 
    hu_raster_files = ["utqa_1t_model_grid.bil", "utqa_2b_model_grid.bil", "utaf_1t_model_grid.bil", "utaf_2b_model_grid.bil", "lta_1t_model_grid.bil", "lta_2b_model_grid.bil"] # ["qa_1t_mb", "qa_2b_mb", "utqa_1t_mb", "utqa_2b_mb"]# #["utqa_1t_mb", "utqa_2b_mb", "utaf_1t_mb", "utaf_2b_mb", "lta_1t_mb", "lta_2b_mb"] #
    #      "utb_1t_model_grid.bil", "utb_2b_model_grid.bil", 
    #      "qa_1t_model_grid.bil", "qa_2b_model_grid.bil",   , "bse_1t_model_grid.bil", "bse_2b.tif_model_grid.bil" 

    #hu_layers = {"qa":["qa_1t_mb", "qa_2b_mb"], "utb":["utb_1t_mb", "utb_2b_mb"], "utqa":["utqa_1t_mb", "utqa_2b_mb"], "utaf":["utaf_1t_mb", "utaf_2b_mb"], "lta":["lta_1t_mb", "lta_2b_mb"], "bse":["bse_1t_mb", "bse_2b_mb.tif"]}
    hu_layers_ordering = ["qa","utb","utqa","utaf","lta","bse"]
    #hu_ext = "_mb"
    vtk_out = 'Campaspe_all_model_mesh'

    out_path = r"C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\model_files\\"

    # Minimum thickness of cells
    min_height = 1
    # Maximum thickness of cells
    max_height = 20

    map_raster_array_to_mesh(hu_raster_path, hu_raster_files, out_path, vtk_out, min_height, max_height)