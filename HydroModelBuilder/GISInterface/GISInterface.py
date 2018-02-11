
class GISInterface(object):

    # Would be good to have something here to clean up old GDAL objects lying around
    # e.g. could have a look through locals() list and

    def __init__(self, GCS=None, PCS=None, VCS=None, pcs_EPSG=None):

        self.geographical_coordinate_system = GCS
        self.projected_coordinate_system = PCS  # Uses osr object from GDAL library
        self.vertical_coordinate_system = VCS
        self.pcs_EPSG = pcs_EPSG

        """
        :param GCS: geographical coordinate system, e.g. GDA_1994_Lambert_Conformal_Conic.
        :param PCS: projected coordinate system, e.g. GDA_1994_MGA_Zone_55.
        :param VCS: vertical coordinate system, e.g. Australian Height Datum (AHD).
        """

    def set_model_boundary_from_polygon_shapefile(self, shapefile_name, shapefile_path=None):
        """
        Function to set the model boundary using a shapefile, assuming polygon type:
        1. For structured mesh this will create a rectangle that creates an envelope on the polygon

        """
        pass

    def set_data_boundary_from_polygon_shapefile(self, shapefile_name, shapefile_path=None, buffer_dist=None):
        """
        Function to set boundary data from which to clip spatial data sets to only include data within the boundary polygon
        """
        pass

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
        pass

    def define_structured_mesh(self, gridHeight, gridWidth):
        """
        Function to define a structured mesh
        """
        pass

    def read_rasters(self, files, path=None):
        """
        Reads in raster files, e.g. for hydrostratigraphy.

        :param files: List of file names for the rasters that are to be read in.
        :param path: Path of the files, which is optional, default path is working directory
        """
        pass

    def map_rasters_to_grid(self, raster_files, raster_path):
        """
        Maps rasters to the defined grid

        :param raster_files: list of raster files to be mapped
        :param raster_path: path of the raster files that are to be mappeed

        writes files to output folder
        """
        pass

    def map_heads_to_mesh_by_layer(self):
        """
        Maps head observations to each layer

        returns array of head valueswith same shape as the mesh
        """
        pass

    def create_basement_bottom(self, hu_raster_path, surface_raster_file, basement_top_raster_file, basement_bot_raster_file, output_path, raster_driver=None):
        """
        Utility to build a bottom basement array where it doesn't exist based on top of bedrock, surface elevation and a thickness function

        writes a raster file for the bottom basement array
        """
        pass

    def build_3D_mesh_from_rasters(self, raster_files, raster_path, minimum_thickness, maximum_thickness):
        """
        Generate the 3D mesh from the
        """
        pass

    def read_polyline(self, filename, path=None):
        """
        Read in polyline data, e.g. for stream definition

        :param filename: filename for the polyline that is to be read in.
        :param path: Path of the files, which is optional, default path is working directory
        """
        pass

    def map_polyline_to_grid(self, polyline_obj):
        """
        Map the polyline object to the grid

        :param polyline_obj
        """
        pass

    def read_points_data(self, filename, path=None):
        """
        Read csv files containing points data, e.g. well locations, screening depths.

        :param filename: filename for the point data that is to be read in.
        :param path: Path of the files, which is optional, default path is working directory
        """
        pass

    def map_points_to_grid(self, points_obj, feature_id=None):
        """
        Read in point shapefile and find which points are in which grid cells

        :param points_obj: point object containing all the shapefile data
        :param model_grid: model grid
        :param feature_id: feature id of interest to be used in returned list, i.e. a list of containing feature id's within each grid cell
        """
        pass

    def map_points_to_3Dmesh(self, points_obj, model_mesh=None, feature_id=None):
        """
        Read in point shapefile and find which points are which mesh cells

        :param points_obj: point object containing all the shapefile data
        :param model_grid: model grid
        :param feature_id: feature id of interest to be used in returned list, i.e. a list of containing feature id's within each mesh cell
        """
        pass

    def getXYpairs(self, points_obj, feature_id=None):
        """
        Return a list of xy pairs for all of the points in the vector layer.

        :param points_obj: points object containing all the shapefile data
        """
        pass
