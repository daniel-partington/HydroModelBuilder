import os
import subprocess
import sys

import numpy as np

# Import the GDAL classes
from osgeo import gdal, gdalconst, ogr, osr

# Import the GIS interface
from HydroModelBuilder.GISInterface.GISInterface import GISInterface

import basement
import create_buffer
# Import spatial processing functions that use GDAL
import fishnet
import map2grid
import map_raster2mesh
import point_values_from_raster
import polygon2points
import raster2polygon
import reproject

class GDALInterface(GISInterface):
    """
    The GDALInterface works within the GISInterface class to apply a range of
    common spatial tasks using the open source GDAL library.

    Functions currently implemented are:
        1. set_model_boundary_from_polygon_shapefile
        2. set_data_boundary_from_polygon_shapefile
        3. set_model_boundary_from_corners
        4. define_structured_mesh
        5. read_rasters
        6. map_rasters_to_grid

         read_polyline
         map_polyline_to_grid

         map_points_to_grid

    Things to fix ...
    inconsistent use of 'mesh' and 'grid', need consistent terminology. Will change
    all to 'mesh' as it is more generic.

    """

    def __init__(self):

        pass

    def _load_shapefile(self, shapefile_name, shapefile_path=None):
        pass

    def set_model_boundary_from_polygon_shapefile(self, shapefile_name, shapefile_path=None):
        """
        Function to set the model boundary using a shapefile, assuming polygon type:
        1. For structured mesh this will create a rectangle that creates an envelope on the polygon

        """
        # Check first that file hasn't been previously processed
        base_name = os.path.join(self.out_data_folder, shapefile_name[:-4])
        new_file = base_name + "_model.shp"

        if os.path.isfile(new_file):
            print 'Using previously generated file: ' + new_file
            driver = ogr.Open(new_file).GetDriver()
            ds = driver.Open(new_file, 0)

            if ds is None:
                print 'Could not open ' + os.path.join(shapefile_path, shapefile_name)
                sys.exit(1)
            # End if
            layer = ds.GetLayer()
            spatialRef = layer.GetSpatialRef()
            # Check spatial ref coords:
            srs = osr.SpatialReference(wkt=spatialRef.ExportToWkt())

        else:
            fname = os.path.join(shapefile_path, shapefile_name)
            driver = ogr.Open(fname).GetDriver()
            ds = driver.Open(fname, 0)

            if ds is None:
                print 'Could not open ' + fname
                sys.exit(1)
            # End if

            layer = ds.GetLayer()
            spatialRef = layer.GetSpatialRef()
            # Check spatial ref coords:
            srs = osr.SpatialReference(wkt=spatialRef.ExportToWkt())

            if self.projected_coordinate_system == None:
                self.projected_coordinate_system = srs

            elif self.projected_coordinate_system == srs:
                print 'No transform required ... continuing'

            else:
                reproject.reproject_layer(ds,
                                          src_cs=srs,
                                          dst_cs=self.projected_coordinate_system,
                                          geom_type=gdal.ogr.wkbPolygon,
                                          copy_dest=new_file)

                ds = driver.Open(base_name + '_model.shp', 0)

                if ds is None:
                    print 'Could not open ' + new_file
                    sys.exit(1)
                # End if

                layer = ds.GetLayer()
            # End if
        # End if

        if self.mesh_type == 'structured':
            for feature in layer:  # should only be 1 feature in the layer for the bounding polygon!!!
                geom = feature.GetGeometryRef()
                xmin, xmax, ymin, ymax = geom.GetEnvelope()
                self.model_boundary = [xmin, xmax, ymin, ymax, self.projected_coordinate_system]
            # End for
        elif self.mesh_type == 'unstructured':
            print 'Mesh type "unstructured" unsupported at the moment'
            sys.exit(1)
        # End if

        ds = None

        self.boundary_poly_file = new_file
        return self.model_boundary, self.boundary_poly_file

    def set_data_boundary_from_polygon_shapefile(self, shapefile_name, shapefile_path=None, buffer_dist=None):

        if shapefile_path == None:
            shapefile_name = shapefile_name
        else:
            shapefile_name = os.path.join(shapefile_path, shapefile_name)
        # End if

        shp_name = os.path.splitext(os.path.basename(shapefile_name))[0]
        new_file = shp_name + "_buffer_" + str(buffer_dist) + ".shp"
        if os.path.isfile(new_file):
            print 'Using previously generated file: ' + new_file
            self.data_boundary = ogr.Open(new_file, 1)
        else:
            self.data_boundary = create_buffer.create_buffer4poly(
                shapefile_name, buffer_distance=buffer_dist)
        # End if
        self.boundary_data_file = new_file
        return self.data_boundary, self.boundary_data_file

    def set_model_boundary_from_corners(self, xmin, xmax, ymin, ymax):
        """
        Function to set model boundary based on x and y bounds: xmin, xmax, ymin, ymax

                    ._________.(xmax, ymax)
                    |         |
                    |         |
                    |         |
                    ._________.
         (xmin,ymin)
        """

        self.model_boundary[0] = xmin
        self.model_boundary[1] = xmax
        self.model_boundary[2] = ymin
        self.model_boundary[3] = ymax

    def define_structured_mesh(self, gridHeight, gridWidth):
        """
        Function to define a structured mesh

        :param gridHeight: Uniform grid height to be used in defining structured mesh
        :param gridWidth: Uniform grid width to be used in defining structured mesh


        """

        new_file = os.path.join(self.out_data_folder,
            'structured_model_grid_{}m'.format(int(gridHeight)), 
            'structured_model_grid_{}m.shp'.format(int(gridHeight)))

        if os.path.isfile(new_file):
            print 'Using previously generated file: ' + new_file
            self.model_mesh = ogr.Open(new_file, 1)
        else:
            structured_mesh = StructuredMesh(xmin=self.model_boundary[0],
                                             xmax=self.model_boundary[1],
                                             ymin=self.model_boundary[2],
                                             ymax=self.model_boundary[3],
                                             gridHeight=gridHeight,
                                             gridWidth=gridWidth)
            self.model_mesh = fishnet.create_fishnet(structured_mesh, self.model_boundary[
                                                     4], copy_dest=self.model_data_folder)
        # End if
        return self.model_mesh

    def read_rasters(self, files, path=None):
        """
        Reads in raster files, e.g. for hydrostratigraphy.

        :param files: List of file names for the rasters that are to be read in.
        :param path: Path of the files, which is optional, default path is working directory
        """
        rasters = {}

#        print path
#        print "============"
#        import sys
#        sys.exit()

        for raster in files:
            rst = os.path.join(self.out_data_folder, raster + '_reproj.bil')
            if os.path.isfile(rst):
                ds = gdal.Open(rst, gdalconst.GA_ReadOnly)
            else:
                ds = gdal.Open(os.path.join(path, raster), gdalconst.GA_ReadOnly)

                if ds is None:
                    print 'Could not open ' + os.path.join(path, raster)
                    sys.exit(1)

                # Check raster GCS
                prj = ds.GetProjection()
                srs = osr.SpatialReference(wkt=prj)
                if self.projected_coordinate_system == srs:
                    print 'No transform required ... continuing'
                else:
                    print 'Reprojecting raster'
                    ds = None

                    target_srs = self.pcs_EPSG  # projected_coordinate_system.ExportToWkt()
                    # FNULL = open(os.devnull, 'w')
                    command = 'gdalwarp -q -t_srs ' + target_srs + ' -r bilinear "' + \
                        path + raster + '" "' + rst + '"'

                    print command  # -dstalpha
                    subprocess.call(command)  # , stdout=FNULL, stderr=subprocess.STDOUT)

                    # ds = reproject.reproject_raster(ds, # raster dataset
                    #                                pixel_spacing=None,
                    #                                src_cs=None,
                    #                                dst_cs=self.projected_coordinate_system,
                    #                                reproj_method=gdal.GRA_Bilinear,
                    #                                create_copy=True,
                    #                                copy_dest=self.out_data_folder + raster + '_reproj.bil',
                    #                                raster_driver="EHdr",
                    #                                set_bounds=None)

                    ds = gdal.Open(rst, gdalconst.GA_ReadOnly)

            print 'Processing: ', raster
            rasters[raster] = [ds.ReadAsArray(), os.path.join(path, raster)]

            ds = None  # close raster file

        return rasters

    def map_rasters_to_grid(self, raster_files, raster_path, raster_driver=None):
        """
        Map the rasters in raster_collection to the model grid

        If the model grid is a structured mesh, then the reproject
        raster can be used to do this mapping

        Interpolation method options in GDAL that can be specified are:
        gdal.GRA_Average
        gdal.GRA_Bilinear
        gdal.GRA_Cubic
        gdal.GRA_CubicSpline
        gdal.GRA_Lanczos
        gdal.GRA_Mode
        gdal.GRA_NearestNeighbour

        """
        simplified_raster_array = {}

        if self.mesh_type == 'structured':
            new_files = [os.path.join(self.out_data_folder_grid, raster +
                                      '_model_grid.bil') for raster in raster_files]
            new_files_exist = [os.path.isfile(os.path.join(self.out_data_folder_grid,
                                                           raster + '_model_grid.bil'))
                               for raster in raster_files]

            if False not in new_files_exist:
                for index, raster in enumerate(new_files):

                    ds = gdal.Open(raster, gdalconst.GA_ReadOnly)  # raster_path

                    if ds is None:
                        print 'Could not open ' + raster
                        sys.exit(1)
                    # End if
                    simplified_raster_array[raster_files[index]] = ds.GetRasterBand(1).ReadAsArray()
                    ds = None
                # end for
            else:
                # layer = self.model_mesh.GetLayer()
                (xmin, xmax, ymin, ymax) = self.model_boundary[0:4]
                # get rows
                rows = int((ymax - ymin) / self.gridHeight) + 1
                # get columns
                cols = int((xmax - xmin) / self.gridHeight) + 1
                pixel_spacing = self.gridHeight
                (xmin, xmax, ymin, ymax) = self.model_mesh.GetLayer().GetExtent()
                set_bounds = (xmin, xmin + cols * self.gridHeight +
                              1.0, ymax - rows * self.gridHeight, ymax)

#                for feature in layer:
#                    geom =feature.GetGeometryRef()
#                    points_string = geom.ExportToWkt()
#                    start = points_string.find('((') + 2
#                    end = points_string.find('))', start)
#                    points_array = [(float(x.split(' ')[0]), float(x.split(' ')[1])) for x in points_string[start:end].split(',')]
#                    pixel_spacing = points_array[1][0] - points_array[0][0]
#                    set_bounds = layer.GetExtent()
#                    break

                # target_srs = self.projected_coordinate_system.ExportToProj4()
                # target_srs = self.projected_coordinate_system.ExportToWkt()
                target_srs = self.pcs_EPSG

                for raster in raster_files:
                    clp_file = os.path.join(self.out_data_folder, raster + '_clipped.bil')
                    # Clip first using boundary polygon
                    # FNULL = open(os.devnull, 'w')
                    command = "gdalwarp -t_srs " + target_srs + ' -cutline "' + self.boundary_poly_file + \
                        '" -crop_to_cutline "' + raster_path + raster + '" "' + clp_file + '"'
                    print command  # -dstalpha
                    subprocess.call(command)  # , stdout=FNULL, stderr=subprocess.STDOUT)

                    # Use clipped raster to map active raster cells to grid
                    ds = gdal.Open(clp_file, gdalconst.GA_ReadOnly)  # raster_path

                    if ds is None:
                        print 'Could not open ' + clp_file
                        sys.exit(1)
                    # End if
                    mapped_raster = reproject.reproject_raster(
                        ds,  # raster_collection[raster], # raster dataset
                        pixel_spacing=pixel_spacing,
                        src_cs=None,
                        dst_cs=self.projected_coordinate_system,
                        reproj_method=gdal.GRA_NearestNeighbour,
                        create_copy=True,
                        copy_dest=os.path.join(self.out_data_folder_grid, raster + '_model_grid.bil'),
                        raster_driver="EHdr",
                                      set_bounds=set_bounds)
                    # for bands in mapped_raster.GetBandCount:
                    simplified_raster_array[raster] = mapped_raster.GetRasterBand(1).ReadAsArray()
                    ds = None
                # End for
            # End if

            return simplified_raster_array

        elif self.mesh_type == 'unstructured':
            pass

        else:
            print 'Something went wrong ... mesh type is not recognised: %s' % self.mesh_type

    def map_heads_to_mesh_by_layer(self):
        """
        Maps head observations to each layer

        returns array of head values with same shape as the mesh
        """
        pass

    def create_basement_bottom(self, hu_raster_path, surface_raster_file, basement_top_raster_file, basement_bot_raster_file, output_path, raster_driver=None):
        """
        Utility to build a bottom basement array where it doesn't exist based on top of bedrock, surface elevation and a thickness function

        writes a raster file for the bottom basement array
        """
        basement.create_basement_bottom(hu_raster_path, surface_raster_file, basement_top_raster_file,
                                        basement_bot_raster_file, output_path, raster_driver=raster_driver)

    def build_3D_mesh_from_rasters(self, raster_files, raster_path, minimum_thickness, maximum_thickness):
        if self.name == None:
            self.name = 'Default'
        return map_raster2mesh.map_raster_array_to_mesh(raster_path, raster_files, self.out_data_folder_grid, self.name + '_model', minimum_thickness, maximum_thickness)

    def read_polyline(self, filename, path=None):
        """
        Read in polyline data, e.g. for stream definition

        :param filename: filename for the polyline that is to be read in.
        :param path: Path of the files, which is optional, default path is working directory
        """
        driver = ogr.GetDriverByName("ESRI Shapefile")
        p_f = os.path.join(path, filename)
        ds = driver.Open(p_f, 0)
        poly_obj = ds.GetLayer()
        if poly_obj == None:
            print 'Could not open ' + p_f
            sys.exit(1)

        srs = poly_obj.GetSpatialRef()
        # print srs.ExportToWkt()
        if srs == self.projected_coordinate_system:
            'No transfrom required ... continue'
        else:
            ds = reproject.reproject_layer(ds,
                                           src_cs=srs,
                                           dst_cs=self.projected_coordinate_system,
                                           create_copy=True,
                                           copy_dest=os.path.join(self.out_data_folder, filename[
                                                                  :-4] + '_model.shp'),
                                           geom_type=ogr.wkbMultiLineString)

        return ds

    def smooth_poly(self, poly_file, smooth_factor=0.0001):
        # Smooth  poly shape using ogr
        # FNULL = open(os.devnull, 'w')
        poly_out = poly_file.split(os.path.sep)[-1]
        poly_out = poly_out.split('.')[0]

        smth = os.path.join(self.out_data_folder, poly_out + '_smoothed.shp')
        smth = smth + '" "'

        # " ".join("ogr2ogr", "-t_srs", smth + poly_file+'"', "-simplify", smooth_factor)
        command = "ogr2ogr -t_srs " + smth + poly_file + '" -simplify ' + smooth_factor
        print command
        subprocess.call(command)  # , stdout=FNULL, stderr=subprocess.STDOUT)

    def map_polyline_to_grid(self, polyline_obj):
        """
        Map the polyline object to the grid

        :param polyline_obj
        :param model_grid

        returns masked array highlighting cells where polylines intersects

        """
        length_and_centroids = map2grid.shp2grid(
            polyline_obj, self.model_mesh, shp_type='poly', data_folder=self.out_data_folder)
        return length_and_centroids

    def polygon2points(self, polygon_obj, to_fname=None, density=10):

        polygon2points.poly2points(polygon_obj, to_fname=to_fname,
                                   density=density, working_directory=self.out_data_folder_grid)
        points_obj = ogr.Open(to_fname)
        points = self.getXYpairs(points_obj)

        return points

    def read_points_data_from_csv(self, filename, path=None):
        """
        Read csv files containing points data, e.g. well locations, screening depths.

        :param filename: filename for the point data that is to be read in.
        :param path: Path of the files, which is optional, default path is working directory

        returns ??
        """

    def read_points_data(self, filename, path=None):
        """
        Read in point data from shapefile, e.g. observation bore locations, rainfall gauges

        :param filename: filename for the point shapefile that is to be read in.
        :param path: Path of the files, which is optional, default path is working directory
        """

        base = os.path.splitext(os.path.basename(filename))[0]
        new_f = base + "_clipped.shp"
        new_file = os.path.join(self.out_data_folder, new_f)

        if os.path.isfile(new_file):
            print 'Using previously generated file: ' + new_file
            driver = ogr.GetDriverByName("ESRI Shapefile")
            ds = driver.Open(new_file, 0)
        else:
            # Clip first using boundary polygon
            target_srs = self.pcs_EPSG

            clipping_poly = self.boundary_data_file

            fn = os.path.join(self.out_data_folder, base) + '_reproj.shp'
            command = 'ogr2ogr -t_srs "' + target_srs + '" "' + fn + '" "' + filename + '"'
            subprocess.call(command)

            command = 'ogr2ogr -clipsrc "' + clipping_poly + '" "' + new_file + '" "' + fn + '"'
            subprocess.call(command)

            driver = ogr.GetDriverByName("ESRI Shapefile")
            ds = driver.Open(new_file, 0)

            """
            if ds == None:
                print 'Could not open ' + new_file
                sys.exit(1)
            # end if

            point_obj = ds.GetLayer()
            if point_obj == None:
                print 'Could not find any layers'
                sys.exit(1)

            srs = point_obj.GetSpatialRef()

            if srs == self.projected_coordinate_system:
                'No transfrom required ... continue'
            else:
                ds = reproject.reproject_layer(ds,
                                               src_cs=srs,
                                               dst_cs=self.projected_coordinate_system,
                                               create_copy=True,
                                               copy_dest=filename[:-4]+'_model.shp', #self.out_data_folder +
                                               geom_type=ogr.wkbMultiLineString)
           # End if
           """
        # End if

        return ds

    def map_points_to_grid(self, points_obj, feature_id=None):
        """
        Map the polyline object to the grid

        :param polyline_obj
        :param model_grid

        returns masked array highlighting cells where polylines intersects

        """
        ids_and_centroids = map2grid.shp2grid(
            points_obj, self.model_mesh, feature_id=feature_id, shp_type='points', data_folder=self.out_data_folder)

        return ids_and_centroids

    def map_points_to_3Dmesh(self, points_obj, model_mesh=None, feature_id=None):
        """
        To be implemented ... don't necessarily need to use GDAL but could be
        based on the k-d tree as in the GWModelBuilder
        """
        pass

    def getXYpairs(self, points_obj, feature_id=1):

        Layer = points_obj.GetLayer()
        Layer.ResetReading()
        # points_list = []
        point_ids = {}
        for feat in Layer:
            geom = feat.GetGeometryRef()
            geom_str = geom.ExportToWkt()
            # points_list += geom_str.split('(')[1].split(')')[0]
            coords = geom_str.split('(')[1].split(')')[0]
            coords = coords.split(' ')
            coords = (float(coords[0]), float(coords[1]))
            point_ids[str(feat.GetField(feature_id))] = coords
        Layer = None
        return point_ids

    def point_values_from_raster(self, points, raster_obj):
        """
        Get values from raster at given points defined by list or by points object
        points must be a list made up of point lists or points shapefile object
        """

        return point_values_from_raster.get_point_values_from_raster(points, raster_obj)

    def map_points_to_raster_layers(self, points, depths, rasters, points_layer):

        p_j = os.path.join
        for index, raster in enumerate(rasters):

            if index % 2 == 1:
                continue

            top = gdal.Open(p_j(self.out_data_folder, rasters[index]))
            bot = gdal.Open(p_j(self.out_data_folder, rasters[index + 1]))

            top_levels = self.point_values_from_raster(points, top)
            bot_levels = self.point_values_from_raster(points, bot)

            # for index, depth in enumerate(depths):
            #    print 't:',top_levels[index],'d:',depth,'b:',bot_levels[index]
            top_levels = np.array(top_levels)
            bot_levels = np.array(bot_levels)
            depths = np.array(depths)
            points_layer[index / 2] = (top_levels > depths) & (bot_levels < depths)

            top = None
            bot = None

        return points_layer
        
    def raster2polygon(self, raster):
        '''
        Function to get the first band of a raster and convert it into 
        a polygon
        '''
        
        return raster2polygon.raster2polygon(raster)
        


class StructuredMesh(object):

    def __init__(self, xmin=None, xmax=None, ymin=None, ymax=None, gridHeight=None, gridWidth=None):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.gridHeight = gridHeight
        self.gridWidth = gridWidth
