import os
import sys
import time

from osgeo import ogr

def poly_line2mesh(poly_line, poly_mesh):

    mesh_layer = poly_mesh.GetLayer()    
    srs = mesh_layer.GetSpatialRef()

    poly_layer = poly_line.GetLayer()    
    
    # Set driver for creating shape files
    driver_inter = ogr.GetDriverByName('ESRI Shapefile')

    polyline_mapped = []

    mesh_layer.ResetReading()

    x = mesh_layer.GetFeatureCount()

    print "Polyline to mesh conversion:"

    for feature in mesh_layer:

        # Track progress because this is very slow
        i = feature.GetFID()        
        sys.stdout.write('\r')
        sys.stdout.write("[%-20s] %d%%" % ('='* int(float(i)/float(x)*20.+1), int(100.0/(x-1.0)*i)))
        sys.stdout.flush()

        grid_cell = feature.GetGeometryRef()
        
        # Setup shape file for each grid cell
        filename = 'mesh_temp.shp'

        while os.path.exists(filename):
            os.remove(filename)

        mesh_layer_ds = None
                
        if os.path.exists('mesh_temp.shp'):
            print 'File still there!!!'
        #end if
        mesh_layer_ds = None
        while mesh_layer_ds == None:
            mesh_layer_ds = driver_inter.CreateDataSource('mesh_temp.shp')

        mesh_temp = None
        while mesh_temp == None:
            mesh_temp = mesh_layer_ds.CreateLayer('mylayer', srs, geom_type=ogr.wkbMultiPolygon)

        featureDefn = mesh_temp.GetLayerDefn()        

        ring = ogr.Geometry(ogr.wkbLinearRing)
        #for i in range(0, grid_cell.GetGeometryCount()):         
        #    print grid_cell.GetPoint(i)            
        (xmin, xmax, ymin, ymax) = grid_cell.GetEnvelope()        
        ring.AddPoint(xmin, ymax)
        ring.AddPoint(xmax, ymax)
        ring.AddPoint(xmax, ymin)
        ring.AddPoint(xmin, ymin)
        ring.AddPoint(xmin, ymax)

        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        # add new geom to layer
        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(poly)
        mesh_temp.CreateFeature(outFeature)
        outFeature.Destroy

        # Setup new shapefile for the intersection
        #fname = 'temp'+str(feature.GetFID())+'.shp'
        fname = 'temp.shp'        
        if os.path.exists(fname):
            os.remove(fname)

        dstshp = None
        while dstshp == None:
            dstshp = driver_inter.CreateDataSource(fname)

        dstlayer = dstshp.CreateLayer('mylayer', srs, geom_type=ogr.wkbLineString)

        poly_layer.Intersection(mesh_temp, dstlayer)
        
        if int(dstlayer.GetFeatureCount()) != 0:        
            length = 0.0
            for feature2 in dstlayer:
                line = feature2.GetGeometryRef()
                length += line.Length()
            #print length, 'm' 
            #print grid_cell.Centroid()        
            centroid_txt = grid_cell.Centroid().ExportToWkt()   
            centroid_nums = centroid_txt.split('(')[1].split(')')[0]           
            centroid = centroid_nums.split(' ')            
            polyline_mapped += [[length, centroid]]           
        
        # clean up
        dstlayer = None            
        dstshp = None
        #if os.path.exists(fname):
        #    os.remove(fname)            
        #driver_inter.DeleteDataSource(fname)
        mesh_temp = None
        mesh_layer_ds = None
        #if os.path.exists('mesh_temp.shp'):
        #    os.remove('mesh_temp.shp')            

    print ""
    
    if os.path.exists('mesh_temp.shp'):
        os.remove('mesh_temp.shp')            
    if os.path.exists(fname):
        os.remove(fname)


    # close shape files
    mesh_layer = None
    poly_layer = None 
    
    return polyline_mapped

if __name__ == "__main__":
    # Open a polyline object
    driver = ogr.GetDriverByName("ESRI Shapefile")        
    ds = driver.Open(r"C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\model_files\Campaspe_Riv_model.shp", 0)
    poly_obj = ds.GetLayer()    
    if poly_obj == None:
        print 'Could not open '
    srs = poly_obj.GetSpatialRef()

    # Open the mesh object
    ds2 = driver.Open(r"C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\model_files\structured_model_grid_20000m\structured_model_grid_20000m.shp", 0)


    mapped_list = poly_line2mesh(ds, ds2)
    for item in mapped_list:
        print item


    ds = None
    ds2 = None