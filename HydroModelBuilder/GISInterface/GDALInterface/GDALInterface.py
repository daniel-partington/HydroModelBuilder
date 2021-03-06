import os
import subprocess
import sys

# Import numpy for numerical and arrays
import numpy as np
# Import the GDAL classes
from osgeo import gdal, gdalconst, ogr, osr

# Import spatial processing functions that use GDAL
import basement
import create_buffer
import fishnet
import map2grid
import map_raster2mesh
import point_values_from_raster
import polygon2points
import polygon2raster
import raster2polygon
import reproject
# Import the GIS interface
from HydroModelBuilder.GISInterface.GISInterface import GISInterface


class GDALInterface(GISInterface):
    """The GDALInterface works within the GISInterface class to apply a range of
    common spatial tasks using the open source GDAL library.

    Things to fix ...
    inconsistent use of 'mesh' and 'grid', need consistent terminology. Will change
    all to 'mesh' as it is more generic."""

    def __init__(self):
        pass

    def _load_shapefile(self, shapefile_name, shapefile_path=None):
        """

        :param shapefile_name:
        :param shapefile_path:  (Default value = None)
        """
        pass
    # End _load_shapefile()

    def _test_osgeo_load(self, obj, file_name):
        """
        Check to see if osgeo object is valid.

        :param obj: osgeo object, if not loaded properly this will be None.
        :param file_name: str, name of file that was loaded
        """
        if obj is None:
            raise IOError('Could not open {}'.format(file_name))
        # End if
    # End _test_osgeo_load()

    def set_model_boundary_from_polygon_shapefile(self, shapefile_name, shapefile_path=None):
        """
        Set the model boundary using a shapefile, assuming polygon type.
        1. For structured mesh this will create a rectangle that creates an envelope on the polygon

        :param shapefile_name: str, name of shapefile to use
        :param shapefile_path: str, path to shapefile

        :returns: tuple[object], (self.model_boundary, self.boundary_poly_file)
        """
        # Check first that file hasn't been previously processed
        base_name = os.path.join(self.out_data_folder, shapefile_name[:-4])
        new_file = base_name + "_model.shp"

        if os.path.isfile(new_file):
            print 'Using previously generated file: ' + new_file
            driver = ogr.Open(new_file).GetDriver()
            ds = driver.Open(new_file, 0)
            self._test_osgeo_load(ds, new_file)

            layer = ds.GetLayer()
            spatialRef = layer.GetSpatialRef()
            # Check spatial ref coords:
            srs = osr.SpatialReference(wkt=spatialRef.ExportToWkt())

        else:
            fname = os.path.join(shapefile_path, shapefile_name)

            try:
                driver = ogr.Open(fname).GetDriver()
            except AttributeError as e:
                print("Failed loading {}".format(fname))
                raise AttributeError(e)
            # End try

            ds = driver.Open(fname, 0)
            self._test_osgeo_load(ds, fname)

            layer = ds.GetLayer()
            spatialRef = layer.GetSpatialRef()
            # Check spatial ref coords:
            srs = osr.SpatialReference(wkt=spatialRef.ExportToWkt())

            if self.projected_coordinate_system is None:
                self.projected_coordinate_system = srs
            elif self.projected_coordinate_system == srs:
                print 'No transform required ... continuing'
            else:
                reproject.reproject_layer(ds,
                                          src_cs=srs,
                                          dst_cs=self.projected_coordinate_system,
                                          geom_type=gdal.ogr.wkbPolygon,
                                          copy_dest=new_file)

                ds = driver.Open(new_file, 0)

                self._test_osgeo_load(ds, new_file)

                layer = ds.GetLayer()
            # End if
        # End if

        if self._mesh_type == 'structured':
            # should only be 1 feature in the layer for the bounding polygon!!!
            for feature in layer:
                geom = feature.GetGeometryRef()
                xmin, xmax, ymin, ymax = geom.GetEnvelope()
                self.model_boundary = [xmin, xmax, ymin, ymax, self.projected_coordinate_system]
            # End for
        elif self._mesh_type == 'unstructured':
            sys.exit('Mesh type "unstructured" unsupported at the moment')
        # End if

        ds = None

        self.boundary_poly_file = new_file
        return self.model_boundary, self.boundary_poly_file

    def set_data_boundary_from_polygon_shapefile(self, shapefile_name, shapefile_path=None,
                                                 buffer_dist=None):
        """
        :param shapefile_name:
        :param shapefile_path:  (Default value = None)
        :param buffer_dist: Default value = None)
        """
        shapefile_name = shapefile_name if not shapefile_path \
            else os.path.join(shapefile_path, shapefile_name)

        shp_name = os.path.splitext(os.path.basename(shapefile_name))[0]
        new_file = os.path.join(shapefile_path, shp_name + "_buffer_" + str(buffer_dist) + ".shp")
        if os.path.isfile(new_file):
            print 'Using previously generated file: ' + new_file
            self.data_boundary = ogr.Open(new_file, 1)
            self._test_osgeo_load(self.data_boundary, new_file)
        else:
            self.data_boundary = create_buffer.create_buffer4poly(
                shapefile_name, buffile=new_file, buffer_distance=buffer_dist)
            self._test_osgeo_load(self.data_boundary, new_file)
        # End if
        self.boundary_data_file = new_file

        return self.data_boundary, self.boundary_data_file
    # End set_data_boundary_from_polygon_shapefile()

    def set_model_boundary_from_corners(self, xmin, xmax, ymin, ymax):
        """Function to set model boundary based on x and y bounds: xmin, xmax, ymin, ymax::

                    ._________.(xmax, ymax)
                    |         |
                    |         |
                    |         |
                    ._________.
                (xmin,ymin)

        :param xmin: float, bottom left x-position
        :param xmax: float, upper right x-position
        :param ymin: float, bottom left y-position
        :param ymax: float, upper right y-position
        """
        self.model_boundary[0] = xmin
        self.model_boundary[1] = xmax
        self.model_boundary[2] = ymin
        self.model_boundary[3] = ymax
    # End set_model_boundary_from_corners

    def define_structured_mesh(self, gridHeight, gridWidth):
        """Function to define a structured mesh

        :param gridHeight: Uniform grid height to be used in defining structured mesh
        :param gridWidth: Uniform grid width to be used in defining structured mesh

        :returns: `model_mesh` property
        """

        new_file = os.path.join(self.out_data_folder,
                                'structured_model_grid_{}m'.format(int(gridHeight)),
                                'structured_model_grid_{}m.shp'.format(int(gridHeight)))

        if os.path.isfile(new_file):
            print 'Using previously generated file: ' + new_file
            self.model_mesh = ogr.Open(new_file, 1)
            self._test_osgeo_load(self.model_mesh, new_file)
        else:
            structured_mesh = StructuredMesh(xmin=self.model_boundary[0],
                                             xmax=self.model_boundary[1],
                                             ymin=self.model_boundary[2],
                                             ymax=self.model_boundary[3],
                                             gridHeight=gridHeight,
                                             gridWidth=gridWidth)

            self.model_mesh = fishnet.create_fishnet(structured_mesh, self.model_boundary[
                                                     4], copy_dest=self.model_data_folder)
            self._test_osgeo_load(self.model_mesh, new_file)
        # End if

        return self.model_mesh
    # End define_structured_mesh()

    def read_rasters(self, files, path=None):
        """Reads in raster files, e.g. for hydrostratigraphy.

        :param files: List of file names for the rasters that are to be read in.
        :param path: Path of the files. (Default value = None)
        :returns: dict[list], raster array data and filename
        """
        rasters = {}

        for raster in files:
            rst = os.path.join(self.out_data_folder, raster + '_reproj.tif')
            if os.path.isfile(rst):
                ds = gdal.Open(rst, gdalconst.GA_ReadOnly)
                self._test_osgeo_load(ds, rst)
            else:
                if path:
                    ds = gdal.Open(os.path.join(path, raster), gdalconst.GA_ReadOnly)
                else:
                    path = ""
                    ds = gdal.Open(raster, gdalconst.GA_ReadOnly)
                # end if

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
                        os.path.join(path, raster) + '" "' + rst + '"'

                    print command

                    try:
                        subprocess.check_output(command, shell=True)
                    except subprocess.CalledProcessError as e:
                        print e

                    ds = gdal.Open(rst, gdalconst.GA_ReadOnly)
                    self._test_osgeo_load(ds, rst)

            print 'Processing: ', raster
            rasters[raster] = [ds.ReadAsArray(), os.path.join(path, raster)]

            ds = None  # close raster file

        return rasters
    # End read_rasters()

    def get_raster_info(self, fname):
        '''
        Open raster file and create dictionary with raster info.

        :param fname: str, filename to get information from.

        :returns: dict, raster array data and metadata
        '''
        gdal.UseExceptions()
        raster = {}
        raster_data = gdal.Open(fname, 0)
        ulx, x_pixel, junk, uly, junk, y_pixel = raster_data.GetGeoTransform()
        raster_band = raster_data.GetRasterBand(1)
        no_data = raster_band.GetNoDataValue()
        raster['array'] = raster_band.ReadAsArray()
        cols = raster_band.XSize
        rows = raster_band.YSize
        raster['metadata'] = {'cols': cols, 'rows': rows, 'nodata': no_data,
                              'ulx': ulx, 'uly': uly, 'pixel_x': x_pixel, 'pixel_y': -y_pixel}
        return raster
    # End get_raster_info()

    def raster_reproject_by_grid(self, infile, outfile, epsg_to=None,
                                 bounds=None, resample_method=None):
        """

        :param infile:
        :param outfile:
        :param epsg_to: (Default value = None)
        :param bounds: (Default value = None)
        :param resample_method: (Default value = None)
        """
        epsg_to = self.pcs_EPSG if not epsg_to else "EPSG:{}".format(epsg_to)

        if not bounds:
            (xmin, xmax, ymin, ymax) = self.model_mesh.GetLayer().GetExtent()
        else:
            (xmin, xmax, ymin, ymax) = bounds
        # End if

        command = "gdalwarp -overwrite -t_srs {} ".format(epsg_to) + \
                  " -te " + " {0} {1} {2} {3} ".format(xmin, ymin, xmax, ymax) + \
                  '"' + infile + '" "' + outfile + \
                  '" -r {}'.format(resample_method)

        try:
            print(subprocess.check_output(command, shell=True))
        except subprocess.CalledProcessError as e:
            print(e)
        # End try
    # End raster_reproject_by_grid()

    def map_rasters_to_grid(self, raster_files, raster_path, raster_driver=None):
        """Map the rasters in raster_collection to the model grid

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

        :param raster_files: list or str, raster file(s) to map to grid
        :param raster_path: str, path to raster files
        :param raster_driver: GDAL driver to use to load rasters (Default value = None)

        :returns: dict or None, raster data or None if `unstructured`
        """
        if type(raster_files) != list:
            # Assume single file given as string:
            raster_files = [raster_files]

        simplified_raster_array = {}
        if self._mesh_type == 'structured':
            new_files = [os.path.join(self.out_data_folder_grid, raster +
                                      '_model_grid.tif') for raster in raster_files]
            new_files_exist = [os.path.isfile(os.path.join(self.out_data_folder_grid,
                                                           raster + '_model_grid.tif'))
                               for raster in raster_files]

            if False not in new_files_exist:
                for index, raster in enumerate(new_files):

                    ds = gdal.Open(raster, gdalconst.GA_ReadOnly)  # raster_path
                    self._test_osgeo_load(ds, raster)

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

                target_srs = self.pcs_EPSG

                for raster in raster_files:
                    clp_file = os.path.join(self.out_data_folder, raster + '_clipped.tif')
                    # Clip first using boundary polygon
                    # FNULL = open(os.devnull, 'w')
                    command = "gdalwarp -overwrite -r bilinear -t_srs " + target_srs + ' -cutline "' + \
                              self.boundary_poly_file + '" -crop_to_cutline "' + \
                              os.path.join(raster_path, raster) + '" "' + clp_file + '"'
                    print command  # -dstalpha

                    try:
                        subprocess.check_output(command, shell=True)
                    except subprocess.CalledProcessError as e:
                        sys.exit(e)

                    # Use clipped raster to map active raster cells to grid
                    ds = gdal.Open(clp_file, gdalconst.GA_ReadOnly)  # raster_path
                    self._test_osgeo_load(ds, clp_file)

                    mapped_raster = reproject.reproject_raster(
                        ds,  # raster_collection[raster], # raster dataset
                        pixel_spacing=pixel_spacing,
                        src_cs=None,
                        dst_cs=self.projected_coordinate_system,
                        reproj_method=gdal.GRA_NearestNeighbour,  # Changed fomr  GRA_Bilinear
                        create_copy=True,
                        copy_dest=os.path.join(self.out_data_folder_grid,
                                               raster + '_model_grid.tif'),  # .bil'),
                        raster_driver="GTiff",  # "EHdr",
                        set_bounds=set_bounds)
                    # for bands in mapped_raster.GetBandCount:
                    simplified_raster_array[raster] = mapped_raster.GetRasterBand(1).ReadAsArray()
                    ds = None
                # End for
            # End if

            return simplified_raster_array

        elif self._mesh_type == 'unstructured':
            pass

        else:
            print('Something went wrong ... mesh type is not recognised: {}'.foramt(self._mesh_type))
        # End if

    # End map_rasters_to_grid

    def get_raster_array(self, raster_fname, band=1):
        """
        Load raster data as an array

        :param raster_fname: str, filename of raster
        :param band: int, band id to load (Default value = 1)

        :returns: numpy array, raster data
        """
        ds = gdal.Open(raster_fname)
        raster_band = ds.GetRasterBand(band)
        array = raster_band.ReadAsArray()
        return array
    # End get_raster_array()

    def map_raster_to_regular_grid_return_array(self, raster_fname,
                                                resample_method='near'):
        """
        :param raster_fname:
        :param resample_method:  (Default value = 'near')
        """

        if not os.path.exists(raster_fname):
            raise IOError("File {} does not exist".format(raster_fname))
        # End if

        (xmin, xmax, ymin, ymax) = self.model_mesh.GetLayer().GetExtent()
        epsg_to = self.pcs_EPSG
        pixel_spacing = self.gridHeight
        if os.path.sep in raster_fname:
            fname = raster_fname.split(os.path.sep)[-1]
        else:
            fname = raster_fname
        # End if

        new_raster_fname = os.path.join(self.out_data_folder_grid, fname[:-4] + '_model.tif')
        if os.path.exists(new_raster_fname):
            print("Found previously generated file: {}".format(new_raster_fname))
        else:
            command = "gdalwarp -overwrite -t_srs {} ".format(epsg_to) + \
                      " -tr {} {} ".format(str(pixel_spacing), str(pixel_spacing)) + \
                      " -te " + " {0} {1} {2} {3} ".format(xmin, ymin, xmax, ymax) + \
                '"' + raster_fname + '" "' + new_raster_fname + \
                '" -r {}'.format(resample_method)

            try:
                print(subprocess.check_output(command, shell=True))
            except subprocess.CalledProcessError as e:
                sys.exit(e)
            # End try
        # End if

        return self.get_raster_array(new_raster_fname)
    # End map_raster_to_regular_grid_return_array()

    def map_heads_to_mesh_by_layer(self):
        """Maps head observations to each layer

        :returns: array of head values with same shape as the mesh
        """
        pass
    # End map_heads_to_mesh_by_layer()

    def create_basement_bottom(self, hu_raster_path, surface_raster_file, basement_top_raster_file,
                               basement_bot_raster_file, output_path, raster_driver=None):
        """Utility to build a bottom basement array where it doesn't exist based on top of bedrock,
        surface elevation and a thickness function.

        writes a raster file for the bottom basement array

        :param hu_raster_path:
        :param surface_raster_file:

        :param basement_top_raster_file:
        :param basement_bot_raster_file:

        :param output_path:
        :param raster_driver:  (Default value = None)
        """

        basement.create_basement_bottom(hu_raster_path, surface_raster_file,
                                        basement_top_raster_file, basement_bot_raster_file,
                                        output_path, raster_driver=raster_driver)
    # End create_basement_bottom()

    def build_3D_mesh_from_rasters(self, raster_files, raster_path, minimum_thickness,
                                   maximum_thickness):
        """

        :param raster_files:
        :param raster_path:

        :param minimum_thickness:
        :param maximum_thickness:
        """

        if self.name is None:
            self.name = 'Default'
        return map_raster2mesh.map_raster_array_to_mesh(raster_path, raster_files,
                                                        self.out_data_folder_grid,
                                                        self.name + '_model',
                                                        minimum_thickness,
                                                        maximum_thickness)
    # End build_3D_mesh_from_rasters()

    def read_poly(self, filename, path="", poly_type='polyline', new_filename_only=False):
        """Read in shapefile data, e.g. for stream definition

        :param filename: filename for the polyline that is to be read in.
        :param path: Path of the files. (Default value = "")
        :param poly_type: str, name of polygon type. (Default value = 'polyline')
        :param new_filename_only: bool, return new filename or the data (Default value = False)
        """
        if poly_type == 'polyline':
            geom_type = ogr.wkbMultiLineString
        if poly_type == 'polygon':
            geom_type = ogr.wkbMultiPolygon
        # End if

        if os.path.sep in filename:
            filename_only = filename.split(os.path.sep)[-1]
        else:
            filename_only = filename
        # End if

        driver = ogr.GetDriverByName("ESRI Shapefile")
        p_f = os.path.join(path, filename)

        ds = driver.Open(p_f, 0)
        self._test_osgeo_load(ds, p_f)

        poly_obj = ds.GetLayer()

        srs = poly_obj.GetSpatialRef()
        if srs == self.projected_coordinate_system:
            print('No transfrom required ... continue')
        else:
            copy_dest = os.path.join(self.out_data_folder, filename_only[:-4] + '_model.shp')
            ds = reproject.reproject_layer(ds,
                                           src_cs=srs,
                                           dst_cs=self.projected_coordinate_system,
                                           create_copy=True,
                                           copy_dest=copy_dest,
                                           geom_type=geom_type)
            self._test_osgeo_load(ds, copy_dest)
        # End if

        ds_copy = ogr.GetDriverByName("Memory").CopyDataSource(
            ds, ds.GetDescription())

        if srs != self.projected_coordinate_system:
            if new_filename_only:
                return copy_dest
            else:
                return ds_copy
            # End if
        else:
            if new_filename_only:
                return os.path.join(path, filename)
            else:
                return ds_copy
            # End if
        # End if
    # End read_poly()

    def smooth_poly(self, poly_file, smooth_factor=0.0001):
        """
        :param poly_file:
        :param smooth_factor:  (Default value = 0.0001)
        """
        # Smooth  poly shape using ogr
        poly_out = poly_file.split(os.path.sep)[-1]
        poly_out = poly_out.split('.')[0]
        smth = os.path.join(self.out_data_folder, poly_out + '_smoothed.shp')
        smth = smth + '" "'
        command = "ogr2ogr -t_srs " + smth + poly_file + '" -simplify ' + smooth_factor
        print(command)
        try:
            subprocess.check_output(command, shell=True)  # , stdout=FNULL, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print(e)
        # End try
    # End smooth_poly()

    def map_polyline_to_grid(self, polyline_obj):
        """Map the polyline object to the grid

        :param polyline_obj
        :param model_grid

        :returns: masked array, highlighting cells where polylines intersects
        """

        length_and_centroids = map2grid.shp2grid(
            polyline_obj, self.model_mesh, shp_type='poly', data_folder=self.out_data_folder)
        return length_and_centroids
    # End map_polyline_to_grid()

    def map_polygon_to_grid(self, polygon_obj, out_fname, pixel_size=None,
                            bounds=None, feature_name=None, field_type=ogr.OFTInteger):
        """Map the polygon object to the grid

        :param polygon_obj:
        :param out_fname:
        :param pixel_size: (Default value = None)
        :param bounds: (Default value = None)
        :param feature_name: (Default value = None)
        :param field_type: (Default value = ogr.OFTInteger)

        :returns: array and dict. Array contings integers corresponding to different features from polygons.
                  The dict contains map of array integers
        """

        return polygon2raster.array_from_rasterize(polygon_obj, out_fname=out_fname,
                                                   pixel_size=pixel_size, bounds=bounds,
                                                   feature_name=feature_name, field_type=field_type)
    # End map_polygon_to_grid()

    def polygon2points(self, polygon_obj, to_fname=None, density=10):
        """

        :param polygon_obj:
        :param to_fname:  (Default value = None)
        :param density: (Default value = 10)

        :returns:
        """

        polygon2points.poly2points(polygon_obj, to_fname=to_fname,
                                   density=density, working_directory=self.out_data_folder_grid)
        points_obj = ogr.Open(to_fname)
        points = self.getXYpairs(points_obj)

        return points
    # End polygon2points()

    def read_points_data_from_csv(self, filename, path=None):
        """Read csv files containing points data, e.g. well locations, screening depths.

        :param filename: filename for the point data that is to be read in.
        :param path: Path of the file. (Default value = None)
        """
        raise NotImplementedError("Method not implemented yet")
    # End read_points_data_from_csv()

    def read_points_data(self, filename, path=None):
        """Read in point data from shapefile, e.g. observation bore locations, rainfall gauges

        :param filename: filename for the point shapefile that is to be read in.
        :param path: Path of the files. (Default value = None)

        :returns: GDAL file pointer
        """
        base = os.path.splitext(os.path.basename(filename))[0]
        new_f = base + "_clipped.shp"
        new_file = os.path.join(self.out_data_folder, new_f)
        driver = ogr.GetDriverByName("ESRI Shapefile")

        if os.path.isfile(new_file):
            print 'Using previously generated file: ' + new_file
        else:
            # Clip first using boundary polygon
            target_srs = self.pcs_EPSG

            clipping_poly = self.boundary_data_file

            fn = os.path.join(self.out_data_folder, base + '_reproj.shp')
            command = 'ogr2ogr -t_srs "' + target_srs + '" "' + fn + '" "' + filename + '"'
            try:
                print(subprocess.check_output(command, shell=True))
            except subprocess.CalledProcessError as e:
                print("stdout output on error:\n" + e.output)

            command = 'ogr2ogr -clipsrc "' + clipping_poly + '" "' + new_file + '" "' + fn + '" -f "ESRI Shapefile"'
            try:
                print(subprocess.check_output(command, shell=True))
            except subprocess.CalledProcessError as e:
                print("stdout output on error:\n" + e.output)

        # End if

        ds = driver.Open(new_file, 0)
        self._test_osgeo_load(ds, new_file)
        return ds
    # End read_points_data()

    def map_points_to_grid(self, points_obj, feature_id=None):
        """Map the points object to the grid

        :param points_obj:
        :param feature_id:  (Default value = None)

        :returns: masked array highlighting cells where polylines intersects
        """

        ids_and_centroids = map2grid.shp2grid(points_obj, self.model_mesh, feature_id=feature_id,
                                              shp_type='points', data_folder=self.out_data_folder)

        return ids_and_centroids

    def map_points_to_3Dmesh(self, points_obj, model_mesh=None, feature_id=None):
        """To be implemented ... don't necessarily need to use GDAL but could be
        based on the k-d tree as in the GWModelBuilder.

        :param points_obj:
        :param model_mesh:  (Default value = None)
        :param feature_id: (Default value = None)
        """

        raise NotImplementedError("Method not implemented yet")
    # End map_points_to_3Dmesh()

    def getXYpairs(self, points_obj, feature_id=1):
        """
        :param points_obj:
        :param feature_id:  (Default value = 1)
        """
        Layer = points_obj.GetLayer()
        Layer.ResetReading()
        point_ids = {}
        for feat in Layer:
            geom = feat.GetGeometryRef()
            geom_str = geom.ExportToWkt()
            coords = geom_str.split('(')[1].split(')')[0]
            coords = coords.split(' ')
            coords = (float(coords[0]), float(coords[1]))
            point_ids[str(feat.GetField(feature_id))] = coords
        Layer = None
        return point_ids
    # End getXYpairs()

    def point_values_from_raster(self, points, raster_obj):
        """Get values from raster at given points defined by list or by points object
        points must be a list made up of point lists or points shapefile object.

        :param points:
        :param raster_obj:
        """
        return point_values_from_raster.get_point_values_from_raster(points, raster_obj)
    # End point_values_from_raster()

    def create_line_string_from_points(self, points, file_in, file_out,
                                       driver_name='ESRI Shapefile'):
        """Convert set of points into a linestring and write to shapefile

        :param points:
        :param file_in:
        :param file_out:
        :param driver_name:  (Default value = 'ESRI Shapefile')
        """

        shp1 = ogr.Open(file_in)
        layer1 = shp1.GetLayer(0)
        srs = layer1.GetSpatialRef()

        driver = ogr.GetDriverByName(driver_name)
        dstshp = driver.CreateDataSource(file_out)
        dstlayer = dstshp.CreateLayer('mylayer', srs, geom_type=ogr.wkbLineString)

        line = ogr.Geometry(ogr.wkbLineString)
        for point in points:
            line.AddPoint(*point)
        # End for

        featureDefn = dstlayer.GetLayerDefn()
        lineFeature = ogr.Feature(featureDefn)
        lineFeature.SetGeometry(line)
        dstlayer.CreateFeature(lineFeature)
        lineFeature.Destroy
        dstshp.Destroy
        dstshp = None
    # End create_line_string_from_points()

    def polyline_explore(self, poly_file):
        '''
        Extract all points from a polyline and that returns a dict containing all points within each feature.

        :param poly_file: str or ogr object, polygon file to use.

        :returns: tuple, (list of geometry points, dict of geometry points)
        '''

        if type(poly_file) == str:
            poly_ds = ogr.Open(poly_file)
        else:
            poly_ds = poly_file

        poly_layer = poly_ds.GetLayer()
        count = poly_layer.GetFeatureCount()
        points_dict = {}
        points_all = []
        point_count = 0
        for index, feat in enumerate(poly_layer):
            geom = feat.GetGeometryRef()
            point_count += geom.GetPointCount()
            if index in range(count):
                points_dict[index] = geom.GetPoints()
                points_all += geom.GetPoints()
            # End if
        # End for

        return points_all, points_dict
    # End polyline_explore()

    def _linestring2points(self, linestring):
        """
        Convert a string description of line points to numerical list.

        :param linestring: str, line type descriptor to convert into list of points

        :returns: list[float], points that describe a line
        """

        linestring = linestring.split("LINESTRING (")[1].split(")")[0].split(',')
        linestring = [[float(x.split(" ")[0]), float(x.split(" ")[1])] for x in linestring]
        return linestring
    # End _linestring2points()

    def polyline_explore_2(self, poly_file):
        """Extract points from multi-linestring

        :param poly_file:
        """

        poly_ds = ogr.Open(poly_file)
        poly_layer = poly_ds.GetLayer()
        points_all = []
        for index, feat in enumerate(poly_layer):
            geom = feat.GetGeometryRef()
            for geo in geom:
                points_all += self._linestring2points(geo.ExportToWkt())

        return points_all
    # End polyline_explore_2()

    def map_points_to_raster_layers(self, points, depths, rasters, points_layer):
        """

        :param points:
        :param depths:

        :param rasters:
        :param points_layer:
        """
        p_j = os.path.join
        for index, raster in enumerate(rasters):

            if index % 2 == 1:
                continue

            top = gdal.Open(p_j(self.out_data_folder, rasters[index]))
            bot = gdal.Open(p_j(self.out_data_folder, rasters[index + 1]))

            top_levels = self.point_values_from_raster(points, top)
            bot_levels = self.point_values_from_raster(points, bot)

            top_levels = np.array(top_levels)
            bot_levels = np.array(bot_levels)
            depths = np.array(depths)
            points_layer[index / 2] = (top_levels > depths) & (bot_levels < depths)

            top = None
            bot = None

        return points_layer
    # End map_points_to_raster_layers()

    def raster2polygon(self, raster):
        '''
        Convert first band of a raster into a polygon.

        :param raster:

        :returns: np.ndarray, raster data
        '''
        return raster2polygon.raster2polygon(raster)
    # End raster2polygon()

# End GDALInterface()


class StructuredMesh(object):
    """TODO: Docs"""

    def __init__(self, xmin=None, xmax=None, ymin=None, ymax=None, gridHeight=None, gridWidth=None):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.gridHeight = gridHeight
        self.gridWidth = gridWidth
