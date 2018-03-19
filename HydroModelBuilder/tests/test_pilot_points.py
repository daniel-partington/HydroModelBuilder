from HydroModelBuilder.Utilities.PilotPoints.pilotpoints import *
import numpy as np
import os
from HydroModelBuilder.GWModelManager import GWModelManager

resolution = 1000
zone_map = {1: 'qa', 2: 'utb', 3: 'utqa', 4: 'utam', 5: 'utaf', 6: 'lta', 
            7: 'bse'}
HGU_map = {'bse':'Bedrock', 'utb':'Newer Volcanics Basalts', 
           'utaf':'Calivil', 'lta':'Renmark', 
           'qa':'Coonambidgal', 'utqa':'Shepparton',
           'utam':'Loxton-Parilla Sands'}


# Get an example mesh from a previous model build
MM = GWModelManager()
model_folder = r'C:/Workspace/part0075/MDB modelling/testbox/00_Campaspe_Cascade/01_steady_state/structured_model_grid_{}m'.format(resolution)
MM.load_GW_model(os.path.join(model_folder, r"01_steady_state_packaged.pkl"))
name = MM.GW_build.keys()[0]
m = MM.GW_build[name]
mesh_array = m.model_mesh3D
zones = len(np.unique(mesh_array[1])) - 1
cell_centers = m.model_mesh_centroids
model_boundary = m.model_boundary

# Generate pilot points from existing mesh with zones and cell centers
# allow for assignment of pilot points at every nth active zone cell using
# skip function.

pp = PilotPoints(output_directory=model_folder)

if resolution == 1000:
    skip=[3, 0, 6, 0, 6, 6, 6] 
    #skip_active=[49, 20, 0, 34, 0, 0, 0]
    skip_active=[3, 20, 0, 34, 0, 0, 0]
elif resolution == 500:
    skip=[0, 0, 12, 0, 12, 12, 12] 
    skip_active=[100, 40, 0, 70, 0, 0, 0]
elif resolution == 20000:
    skip=[0, 0, 0, 0, 0, 0, 0] 
    skip_active=[0, 0, 0, 0, 0, 0, 0]
else:
    skip=[0,  0, 3, 0, 2, 3, 3] 
    skip_active=[3, 20, 0, 4, 0, 0, 0]



a = pp.generate_points_from_mesh(mesh_array, cell_centers, 
               skip=skip, 
        skip_active=skip_active,
        #zone_prop_dict={0:30.0, 1:2.0, 2:2.0, 3:50.0, 4:45.0, 5:30.0, 6:0.01},
        zone_prop_dict={0:0.01, 1:0.02, 2:0.4, 3:0.5, 4:0.6, 5:0.2, 6:0.10},
        add_noise=True
        )

pp.write_settings_fig()
pp.write_grid_spec(mesh_array, model_boundary, delc=resolution, delr=resolution)
pp.write_struct_file(mesh_array, nugget=0.0, 
                  transform='log',numvariogram=1, variogram=0.15, 
                  vartype=1, bearing=0.0, a=10000.0, anisotropy=1.0)

if resolution == 1000:
    search_radius = [40000, 20000, 15000, 20000, 20000, 20000, 20000]
else:
    search_radius = [30000, 20000, 40000, 20000, 40000, 50000, 20000]

pp.setup_pilot_points_by_zones(mesh_array, zones, search_radius)    

pp.generate_cov_mat_by_zones(zones)

pp.run_pyfac2real_by_zones(zones)



hk = pp.val_array

nrow, ncol = mesh_array[1][0].shape[0], mesh_array[1][0].shape[1]    
delc, delr = [resolution]*2
x0, y0 = model_boundary[0], model_boundary[3]
grid_info = [ncol, nrow, delc, delr, x0, y0]
#mesh3DToVtk(mesh_array, grid_info, hk, 'hk', '', 'hk')
    
###########################################################################
from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface
data_folder = r"C:/Workspace/part0075/MDB modelling/testbox/PEST5000/master"
modflow_model = flopyInterface.ModflowModel(MM.GW_build[name], data_folder=data_folder)
modflow_model.buildMODFLOW()
import flopy
import matplotlib.pyplot as plt
import os
import time

params = False
if params:
       
    default_array = mesh_array[1].astype(float)
    Zone = np.copy(default_array)
    Kh = np.copy(default_array)
    Kv = np.copy(default_array)
    Sy = np.copy(default_array)
    SS = np.copy(default_array)
    
    def create_pp_points_dict(zone_map, Zone, prop_array, prop_name, m):
        points_values_dict = {}
        for index, key in enumerate(zone_map.keys()):
            for index2, param in enumerate(m.parameters.param_set[prop_name + zone_map[key]]):
                if index2 == 0:
                    points_values_dict[index] = [m.parameters.param[param]['PARVAL1']]
                else: 
                    points_values_dict[index] += [m.parameters.param[param]['PARVAL1']]
        return points_values_dict    
        
    def update_pilot_points(zone_map, Zone, prop_array, par_name, prop_name, prop_folder, m, prop_array_fname):
        points_values_dict = create_pp_points_dict(zone_map, Zone, prop_array, prop_name, m)
        p = m.pilot_points[par_name]
        zones = len(zone_map.keys())
        p.output_directory = os.path.join(data_folder, prop_folder)
        p.update_pilot_points_files_by_zones(zones, points_values_dict)
        time.sleep(3)
        p.run_pyfac2real_by_zones(zones) 
        #p.save_mesh3D_array(filename=os.path.join(data_folder, prop_array_fname))
        return p.val_array
    
    Kh = update_pilot_points(zone_map, Zone, Kh, 'hk', 'kh', 'hk_pilot_points',
                             m, 'hk_val_array')  
    Kv = Kh * 0.1

    Sy = update_pilot_points(zone_map, Zone, Sy, 'sy', 'sy', 'sy_pilot_points',
                             m, 'sy_val_array')
    
    SS = update_pilot_points(zone_map, Zone, SS, 'ss', 'ss', 'ss_pilot_points',
                             m, 'ss_val_array')





b = a[0]
# First step is to set up the plot
width = 20
height = 10
multiplier = 1.
fig = plt.figure(figsize=(width * multiplier, height * multiplier))

for key in b.keys():
    fname = os.path.join(pp.output_directory, 'values{}.ref'.format(key))
    with open(fname, 'r') as f:
        lines = f.readlines()
        vals = []
        for index, line in enumerate(lines):
            #if index == 0: continue
            temp = line.strip('\n').split()
            vals += [float(x) for x in line.strip('\n').split()]
    
    #vmin = 1
    #vmax = 2                
    vals = np.array(vals)
    vals = np.reshape(vals, mesh_array[0][0].shape)
    #fig = plt.figure()
    ax = fig.add_subplot(2, 4, key + 1, aspect='equal')
    ax.set_title(HGU_map[zone_map[key + 1]])
    modelmap = flopy.plot.ModelMap(model=modflow_model.mf) 
    array = modelmap.plot_array(a[3][key], masked_values=[0], alpha=0.8, cmap='gray', vmin=0, vmax=2)
    #modelmap.plot_grid()    
    #array = modelmap.plot_array(hk[key], masked_values=[1E6, 0], alpha=1.0)#, 
    array = modelmap.plot_array(hk[key], masked_values=[0], alpha=1.0, cmap='viridis')#, 
                                #vmin=vmin, vmax=vmax)
    #array = modelmap.plot_array(vals, masked_values=[1E6], alpha=1.0)#, 
                                #vmin=vmin, vmax=vmax)
    start, end = ax.get_xlim()
    start = start // 1000 * 1000 + 1000
    end = end // 1000 * 1000 - 1000
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    plt.scatter([x[0] for x in b[key]], [x[1] for x in b[key]], facecolors='none', color='black', alpha=0.5)
    cbar_ax = fig.add_axes() #[0.67, 0.055, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax)        

fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.12)

fig = plt.figure(figsize=(width * multiplier, height * multiplier))

for key in b.keys():
    fname = os.path.join(pp.output_directory, 'values{}.ref'.format(key))
    with open(fname, 'r') as f:
        lines = f.readlines()
        vals = []
        for index, line in enumerate(lines):
            #if index == 0: continue
            temp = line.strip('\n').split()
            vals += [float(x) for x in line.strip('\n').split()]
    
    #vmin = 1
    #vmax = 2                
    vals = np.array(vals)
    vals = np.reshape(vals, mesh_array[0][0].shape)
    #fig = plt.figure()
    ax = fig.add_subplot(2, 4, key + 1, aspect='equal')
    ax.set_title(HGU_map[zone_map[key + 1]])
    modelmap = flopy.plot.ModelMap(model=modflow_model.mf) 
    #array = modelmap.plot_array(a[3][key], masked_values=[0], alpha=0.8, cmap='gray', vmin=0, vmax=2)
    #modelmap.plot_grid()    
    #array = modelmap.plot_array(hk[key], masked_values=[1E6, 0], alpha=1.0)#, 
    current = np.ma.masked_where(~a[3][key], hk[key])
    print("maximum: {}, minimum: {}".format(np.amax(current), np.amin(current)))
    array = modelmap.plot_array(np.ma.masked_where(~a[3][key], hk[key]), masked_values=[0], alpha=1.0, cmap='viridis')#, 
                                #vmin=vmin, vmax=vmax)
    #array = modelmap.plot_array(vals, masked_values=[1E6], alpha=1.0)#, 
                                #vmin=vmin, vmax=vmax)
    start, end = ax.get_xlim()
    start = start // 1000 * 1000 + 1000
    end = end // 1000 * 1000 - 1000
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    plt.scatter([x[0] for x in b[key]], [x[1] for x in b[key]], facecolors='none', color='black', alpha=0.5)
    cbar_ax = fig.add_axes() #[0.67, 0.055, 0.01, 0.42])
    fig.colorbar(array, cax=cbar_ax)        

fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.12)


from matplotlib import colors

cmap = colors.ListedColormap(['lightgrey', 'white'])

fig = plt.figure(figsize=(width * multiplier, height * multiplier))

total_pp = 0
for key in b.keys():
    print("{}:{}".format(key, len(b[key])))
    total_pp += len(b[key])
print total_pp

counter = 0
for key in b.keys():
    if key in [1, 3, 6]:
        continue
    fname = os.path.join(pp.output_directory, 'values{}.ref'.format(key))
    with open(fname, 'r') as f:
        lines = f.readlines()
        vals = []
        for index, line in enumerate(lines):
            #if index == 0: continue
            temp = line.strip('\n').split()
            vals += [float(x) for x in line.strip('\n').split()]
    
    #vmin = 1
    #vmax = 2                
    vals = np.array(vals)
    vals = np.reshape(vals, mesh_array[0][0].shape)
    #fig = plt.figure()
    ax = fig.add_subplot(1, 4, counter + 1, aspect='equal')
    ax.set_title(HGU_map[zone_map[key + 1]])
    modelmap = flopy.plot.ModelMap(model=modflow_model.mf) 
    #array = modelmap.plot_array(a[3][key], masked_values=[0], alpha=0.8, cmap='gray', vmin=0, vmax=2)
    #modelmap.plot_grid()    
    #array = modelmap.plot_array(hk[key], masked_values=[1E6, 0], alpha=1.0)#, 
    array = modelmap.plot_array(~a[3][key], masked_values=[], alpha=1.0, cmap=cmap)#, 
                                #vmin=vmin, vmax=vmax)
    #array = modelmap.plot_array(vals, masked_values=[1E6], alpha=1.0)#, 
                                #vmin=vmin, vmax=vmax)
    start, end = ax.get_xlim()
    start = start // 1000 * 1000 + 1000
    end = end // 1000 * 1000 - 1000
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    if counter > 0:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel('Northing')
    if counter == 1:
        ax.set_xlabel('                                                 Easting', labelpad=20)
       
    plt.scatter([x[0] for x in b[key]], [x[1] for x in b[key]], facecolors='black', color='black', alpha=1.0)
    #cbar_ax = fig.add_axes() #[0.67, 0.055, 0.01, 0.42])
    #fig.colorbar(array, cax=cbar_ax)        
    counter += 1
fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.12)