import os
import time
import sys

from osgeo import ogr


def points2mesh(points, poly_mesh, feature_id=None):

    mesh_layer = poly_mesh.GetLayer()    
    srs = mesh_layer.GetSpatialRef()

    points_layer = points.GetLayer()    
    
    # Set driver for creating shape files
    driver_inter = ogr.GetDriverByName('ESRI Shapefile')

    points_mapped = []

    mesh_layer.ResetReading()

    x = mesh_layer.GetFeatureCount()

    print "Points to mesh conversion:"

    for feature in mesh_layer:

        # Track progress because this is very slow
        i = feature.GetFID()        
        sys.stdout.write('\r')
        sys.stdout.write("[%-20s] %d%%" % ('='* int(float(i)/float(x)*20.+1), int(100.0/(x-1.0)*i)))
        sys.stdout.flush()
        
        grid_cell = feature.GetGeometryRef()
        
        # Setup shape file for each grid cell
        if os.path.exists('mesh_temp.shp'):
            time.sleep(0.005)
            try:
                os.remove('mesh_temp.shp')
            except:
                print 'Could not delete ... trying again'                
                time.sleep(0.1)
                try: 
                    os.remove('mesh_temp.shp')
                except:
                    print 'Still could not delete ... trying again'                
                    time.sleep(180)
                    os.remove('mesh_temp.shp')

        if os.path.exists('mesh_temp.shp'):
            print 'File still there!!!'
        #end if
        mesh_layer_ds = driver_inter.CreateDataSource('mesh_temp.shp')
        if mesh_layer_ds == None:
            if os.path.exists('mesh_temp.shp'):
                os.remove('mesh_temp.shp')
                mesh_layer_ds = driver_inter.CreateDataSource('mesh_temp.shp')
            
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
        dstshp = driver_inter.CreateDataSource(fname)
        dstlayer = dstshp.CreateLayer('mylayer', srs, geom_type=ogr.wkbPoint)

        points_layer.Intersection(mesh_temp, dstlayer)
        
        if int(dstlayer.GetFeatureCount()) != 0:        
            point_ids = []
            for feature2 in dstlayer:
                point_ids += [feature2.GetField(feature_id)]
            # end for                
            centroid_txt = grid_cell.Centroid().ExportToWkt()   
            centroid_nums = centroid_txt.split('(')[1].split(')')[0]           
            centroid = centroid_nums.split(' ')            
            points_mapped += [[point_ids, centroid]]           
        #end if
            
        # clean up
        dstlayer = None            
        dstshp = None
        mesh_temp = None
        mesh_layer_ds = None

    
    print ""    
    # Clean up
    if os.path.exists('mesh_temp.shp'):
        os.remove('mesh_temp.shp')
    if os.path.exists(fname):
        os.remove(fname)

    # close shape files
    mesh_layer = None
    points_layer = None 
    
    return points_mapped

if __name__ == "__main__":
    # Open a points object
    driver = ogr.GetDriverByName("ESRI Shapefile")        
    ds = driver.Open(r"C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\model_files\pumping wells_clipped.shp", 0)
    poly_obj = ds.GetLayer()    
    if poly_obj == None:
        print 'Could not open '
    srs = poly_obj.GetSpatialRef()

    # Open the mesh object
    ds2 = driver.Open(r"C:\Workspace\part0075\MDB modelling\integrated\Modules\Groundwater\model_files\structured_model_grid_20000m\structured_model_grid_20000m.shp", 0)


    mapped_list = points2mesh(ds, ds2, feature_id = "OLD ID")
    #for item in mapped_list:
    #    print item
    print mapped_list[0]

    ds = None
    ds2 = None