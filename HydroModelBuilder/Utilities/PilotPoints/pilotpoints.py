'''
Module to create files required for creating pilot points with PEST utilities
PPK2FAC and PPKREG and FAC2REAL

Works with mesh object from GWModelBuilder that contains number of rows and
columns as well as the zone for each cell.

Contains functions for generating files:
    structure file
    grid specification
    pilot points
    zonal integer array

'''
import os
import subprocess
import time

import numpy as np

class PilotPoints(object):

    def __init__(self, output_directory=None):
        self.min_allowable_points_separation = 0.0
        self.grid_spec_fname = 'grid.spc'
        self.points_fname = 'points.pts'
        self.zone_fname = 'zone.inf'
        self.struct_fname = 'struct.dat'
        self.ppk2fac_exe = os.path.join(os.getcwd(), 'ppk2fac.exe')
        self.fac2real_exe = 'fac2real.exe'
        self.output_directory = os.path.join(output_directory, 'pilot_points')                                        
        if not os.path.exists(self.output_directory):
            os.mkdir(self.output_directory)
            
        self.num_ppoints_by_zone = {}
        
    def write_settings_fig(self, colrow='yes', date=r'mm/dd/yyyy'):
        '''
        Function to write the "settings.fig" file that is required to be present 
        when running ppk2fac
        '''
        out_fname = os.path.join(self.output_directory, 'settings.fig')
        with open(out_fname, 'w') as f:
            f.write('colrow={} \n'.format(colrow))
            f.write('date={}'.format(date))
    
    @staticmethod
    def _zone_array2layers(zone_array):
        '''
        Function to generate masked layer arrays for each zone
        '''
        zones = [x + 1 for x in range(len(np.unique(zone_array)) - 1)]
        layers = zone_array.shape[0]
        zone_mask2D = {}
        for index, zone in enumerate(zones):
            zone_mask = np.ma.masked_array(zone_array, 
                                           zone_array == zone).mask
        
            zone_mask2D[index] = np.full_like(zone_mask[0], False, dtype=bool)
            for layer in range(layers):
               zone_mask2D[index] |= zone_mask[layer] 
    
        return zone_mask2D        
            
    def generate_points_from_mesh(self, mesh_array, cell_centers, zone_prop_dict=None, 
                                  skip=0, skip_active=0, add_noise=False):
        '''
        Function to generate pilot points based in the zonal array in the mesh 
        array object.
        
        Returns points array, points zone array and initial values.
        '''
    
        # First get the number of zones and layers from the mesh_array object
        zones = [x + 1 for x in range(len(np.unique(mesh_array[1])) - 1)]
        points_dict = {}
        points_zone_dict = {}
        points_val_dict = {}

        zone_mask2D = self._zone_array2layers(mesh_array[1])
        self.zone_array = mesh_array[1]
        # For each zone find the extent through cycling through the layers and 
        # creating a mask.
    
        for index, zone in enumerate(zones):
            if type(skip) == list:
                skip_n = skip[index]
            else:
                skip_n = skip
            if type(skip_active) == list:
                skip_active_n = skip_active[index]
            else:
                skip_active_n = skip_active
                
                
            points_dict[index] = []
            points_zone_dict[index] = []
            points_val_dict[index] = []
    
            row_points = 0
            col_points = 0
            points_active = 0
            for row in range(zone_mask2D[index].shape[0]):
                if skip != 0:
                    row_points += 1
                    if row_points > 1 and row_points < skip_n:
                        continue
                    elif row_points == skip_n:
                        row_points = 0
                        continue
                for col in range(zone_mask2D[index].shape[1]):
                    if skip_n != 0:
                        col_points += 1
                        if col_points > 1 and col_points < skip_n:
                            continue
                        elif col_points == skip_n:
                            col_points = 0
                            continue

                    if zone_mask2D[index][row][col]:
                        if skip_active_n != 0:
                            points_active += 1
                            if points_active > 1 and points_active < skip_active_n:
                                continue
                            elif points_active == skip_active_n:
                                points_active = 0
                                continue
                            
                        points_dict[index] += [(cell_centers[0][row][col], 
                                                cell_centers[1][row][col])]
                        
                        points_zone_dict[index] += [index]
                        if zone_prop_dict is not None:
                            x = zone_prop_dict[index]
                            if add_noise:
                                points_val_dict[index] += [x + 0.1*x*np.random.standard_normal()]
                            else:
                                points_val_dict[index] += [x] 
                        else: 
                            if add_noise:
                                points_val_dict[index] += [1.0 + np.random.lognormal(-0.5, 0.2)]
                            else:    
                                points_val_dict[index] += [1.0]

        self.points_dict = points_dict
        self.points_zone_dict = points_zone_dict
        self.points_val_dict = points_val_dict
        self.zone_mask2D = zone_mask2D
        return points_dict, points_zone_dict, points_val_dict, zone_mask2D
    
    def write_grid_spec(self, mesh_array, model_boundary, 
                        grid_spec_fname='grid.spc', delc=None, delr=None):

        '''
        Function to writes the grid specification file using a mesh array
        for use in the PPK2FAC utility
        
        This only works for model grids that are evenly spaced along the rows 
        and columns
        
        ** Future implementation should allow for any type of array which will need
           to be inherited from the model_mesh object in GWModelBuilder
        '''
        
        layer = 0
        row, col = mesh_array[1][layer].shape[0], mesh_array[1][layer].shape[1]
        with open(os.path.join(self.output_directory, grid_spec_fname), 'w') as f:
            f.write('{0} {1}\n'.format(row, col))
            f.write('{0} {1} {2} \n'.format(model_boundary[0], model_boundary[3], 0.0))
            f.write(" ".join(['{}'.format(delc) for x in range(col)]) + '\n')
            f.write(" ".join(['{}'.format(delr) for x in range(row)]) + '\n')
    
    def write_pilot_points_file(self, pilot_points, pp_zones, initial_values,
                                points_fname='points.pts'):
        '''
        Function to write the pilot point files containing name, easting, 
        northing, zone and value for each pilot point
        '''
        with open(os.path.join(self.output_directory, points_fname), 'w') as f:
            for index, point in enumerate(pilot_points):
                f.write('pp{0} {1} {2} {3} {4} \n'.format(index, point[0], 
                        point[1], pp_zones[index] + 1, initial_values[index]))

    def update_pilot_points_file_by_zone(self, new_values, zone,
                                points_fname='points.pts'):
        '''
        Function to write the pilot point files containing name, easting, 
        northing, zone and value for each pilot point
        '''
        with open(os.path.join(self.output_directory, points_fname), 'w') as f:
            for index, point in enumerate(self.points_dict[zone]):
                f.write('pp{0} {1} {2} {3} {4} \n'.format(index, point[0], 
                        point[1], self.points_zone_dict[zone][index] + 1, new_values[index]))
                
                
    def write_zone_file(self, zone_array, zone_fname='zone.inf'):
        '''
        Function to write the zone array file for use in ppk2fac based on an 
        integer array of cell zones.
        '''
        layer = 0
        row, col = zone_array[layer].shape[0], zone_array[layer].shape[1]
        zone_array = zone_array[int(zone_fname[4])]
        zone_array = zone_array.astype(int)
        zone_array[zone_array == -1] = 0
        zone_array[zone_array != 0] = int(zone_fname[4]) + 1
        with open(os.path.join(self.output_directory, zone_fname), 'w') as f:
            f.write('{0} {1}\n'.format(col, row))
            for row in range(zone_array.shape[0]):
                f.write(" ".join([
                        '{}'.format(x) for x in zone_array[row].flatten()]))
                f.write('\n')
    
    def write_struct_file(self, mesh_array, struct_fname='struct.dat', nugget=0.0, 
                          transform='log',numvariogram=1, variogram=0.15, 
                          vartype=2, bearing=0.0, a=500.0, anisotropy=1.0):
        '''
        Function to write the structure file for use in ppk2fac using default 
        values as specified in the arguments to the function, which can then be
        manually modified throught "struct_file_out"
        '''
        #layer = 0
        #zones = len(np.unique(mesh_array[1][layer])) - 1
        zones = len(np.unique(mesh_array[1])) - 1
        
        with open(os.path.join(self.output_directory, struct_fname), 'w') as f:
            for zone in range(zones):
                f.write('STRUCTURE structure{}\n'.format(zone))
                f.write('  NUGGET {}\n'.format(nugget))
                f.write('  TRANSFORM {}\n'.format(transform))
                f.write('  NUMVARIOGRAM {}\n'.format(numvariogram))
                f.write('  VARIOGRAM vario{0} {1}\n'.format(zone, variogram))
                f.write('END STRUCTURE\n')
                f.write('\n')    
            for zone in range(zones):
                f.write('VARIOGRAM vario{}\n'.format(zone))
                f.write('  VARTYPE {}\n'.format(vartype))
                f.write('  BEARING {}\n'.format(bearing))
                f.write('  A {}\n'.format(a))
                f.write('  ANISOTROPY {}\n'.format(anisotropy))
                f.write('END VARIOGRAM\n')
                f.write('\n')
                
    def write_ppk2fac_instruct(self, mesh_array, grid_spec_fname, points_fname, 
                               struct_fname, zone_fname=None, 
                               instruct_fname='ppk2fac.in',
                               min_allowable_points_separation=0.0, kriging='o', 
                               search_radius=10000, min_pilot_points=1, 
                               max_pilot_points=5, factor_fname='factors.dat', 
                               formatted='f', sd_fname='sd.ref', 
                               reg_fname='reg.dat'):
    
        '''
        Write the instructions to be used with executing ppk2fac 
        to allow batch execution
        '''
        zones = len(np.unique(mesh_array[1])) - 1
        with open(os.path.join(self.output_directory, instruct_fname), 'w') as f:
            f.write('{} \n'.format(grid_spec_fname)) 
            f.write('{} \n'.format(points_fname))
            f.write('{} \n'.format(min_allowable_points_separation))
            if zone_fname is not None:
                f.write('{} \n'.format(zone_fname))
            else:
                f.write(' \n')
            #end if    
            f.write('{} \n'.format(struct_fname))
            f.write('\n') # First zone represents inactive cells so skipped
            if zone_fname is not None:
                    f.write('structure{}\n'.format(int(zone_fname[4])))
                    f.write('{}\n'.format(kriging))
                    f.write('{}\n'.format(search_radius))
                    f.write('{}\n'.format(min_pilot_points))
                    f.write('{}\n'.format(max_pilot_points))
            else:
                for zone in range(zones):
                    f.write('structure{}\n'.format(zone))
                    f.write('{}\n'.format(kriging))
                    f.write('{}\n'.format(search_radius))
                    f.write('{}\n'.format(min_pilot_points))
                    f.write('{}\n'.format(max_pilot_points))
                    
            f.write('{}\n'.format(factor_fname))
            f.write('{}\n'.format(formatted))
            f.write('{}\n'.format(sd_fname))
            f.write('{}'.format(reg_fname))
                
    def run_ppk2fac(self, ppk2fac_exe, instruct_fname=None):
        '''
        Function to run ppk2fac using generated input commands from
        write_ppk2fac_instruct
        '''
        
        # Change into location where input files for ppk2fac are located
        cwd = os.getcwd()
        os.chdir(self.output_directory)
        
        # Run ppk2fac
        if instruct_fname is not None:
            subprocess.Popen([ppk2fac_exe, '<', instruct_fname], shell=True)
        else:
            command = ppk2fac_exe
            subprocess.call(command)
        #end if
        
        # Return to current working directory
        os.chdir(cwd)
            
    def write_fac2real_instruct(self, factors_file, points_file, lower_lim, 
                                upper_lim, out_file, 
                                instruct_fname='fac2real.in'):
        '''
        Write the instructions to be used with executing fac2real
        - Assumes formatted file for factors file
        - Assumes use of singular value rather than array
        - Sets value for elements to which no interpolation takes place
        '''
        with open(os.path.join(self.output_directory, instruct_fname), 'w') as f:
            f.write('{} \n'.format(factors_file)) 
            f.write('f \n')
            f.write('{} \n'.format(points_file))
            f.write('s \n')
            f.write('{} \n'.format(lower_lim))
            f.write('s \n')
            f.write('{} \n'.format(upper_lim))
            f.write('{} \n'.format(out_file))
            f.write('1e35 \n')
    
    def run_fac2real(self, fac2real_exe, instruct_fname=None):
        '''
        Function to run fac2real using generated input commands from
        write_fac2real_instruct
        '''

        # Change into location where input files for ppk2fac are located
        cwd = os.getcwd()
        os.chdir(self.output_directory)

        # Run fac2real
        if instruct_fname is not None:
            subprocess.call([fac2real_exe, '<', instruct_fname], shell=True)
        else:
            command = fac2real_exe
            subprocess.call(command)
        #end if
        
        # Return to current working directory
        os.chdir(cwd)
    
    def run_pyfac2real(self, pp_fname='points.dat', factors_fname='factors.dat', 
                       out_fname='values.ref', lower_lim=1.0E-6, upper_lim=1.0E+6):
        '''
        Function to run Python implementation of fac2real using the implementation
        in pyEmu by Jeremy White, see:
        https://github.com/jtwhite79/pyemu/blob/master/pyemu/utils/gw_utils.py
        '''
        
        import pyemu
        
        pyemu.utils.fac2real(str(os.path.join(self.output_directory, pp_fname)), 
                             str(os.path.join(self.output_directory, factors_fname)), 
                             out_file=os.path.join(self.output_directory, out_fname), 
                             upper_lim=upper_lim, lower_lim=lower_lim)

    def _create_mesh3D_array_from_values(self, zone_array, values_dir='', out_dir=''):
        '''
        Function to read in all of the different values arrays for each layer and
        construct the 3D array pertaining to pilot point values as defined by their 
        zones
        '''
        val_array = np.zeros_like(zone_array)
        zones = [x + 1 for x in range(len(np.unique(zone_array)) - 1)]
                 
        for index, zone in enumerate(zones):
            fname = os.path.join(self.output_directory, 'values{}.ref'.format(index))
            with open(fname, 'r') as f:
                lines = f.readlines()
                vals = []
                for ind, line in enumerate(lines):
                    vals += [float(x) for x in line.strip('\n').split()]
            
            vals = np.array(vals)
            vals = np.reshape(vals, zone_array[0].shape)
            for layer in range(zone_array.shape[0]):
                val_array[layer][zone_array[layer]==[zone]] = vals[zone_array[layer]==[zone]]

        self.val_array = val_array           
        #return val_array
        
    def save_mesh3D_array(self, filename='val_array', data_format='binary'):
        '''
        Function to save the created 3D mesh array to file with filename 
        as 'filename' and as either 'data_format' as 'binary' or 'ascii'
        '''
        if data_format == 'ascii':
            np.savetxt(filename, self.val_array)
        elif data_format == 'binary':
            np.save(filename, self.val_array)
        else:
            print 'Data format not recognised, use "binary" or "ascii"'
        # end if

    def setup_pilot_points_by_zones(self, mesh_array, zones, search_radius):
        '''
        Function to set up pilot points files, zone files, ppk2fac instructions,
        and to run ppk2fac to create the factors files required for fac2real.
        
        This is based on creating a pilot point grid per zone in the mesh.
        '''
        for zone in range(zones):
            self.write_pilot_points_file(self.points_dict[zone], 
                                         self.points_zone_dict[zone], 
                                         self.points_val_dict[zone], 
                                         points_fname='points{}.pts'.format(zone))

            self.num_ppoints_by_zone[zone] = len(self.points_dict[zone])            

            self.write_zone_file(self.zone_mask2D, 
                                 zone_fname='zone{}.inf'.format(zone))
    
            self.write_ppk2fac_instruct(mesh_array, 'grid.spc',
                                   'points{}.pts'.format(zone), 
                                   zone_fname='zone{}.inf'.format(zone),
                                   instruct_fname = 'ppk2fac{}.in'.format(zone),
                                   struct_fname='struct.dat', 
                                   factor_fname='factors{}.dat'.format(zone), 
                                   sd_fname='sd{}.ref'.format(zone), 
                                   reg_fname='reg{}.dat'.format(zone),
                                   min_allowable_points_separation=0.0, kriging='o', 
                                   search_radius=search_radius[zone], 
                                   min_pilot_points=1, 
                                   max_pilot_points=50)
            
            #if os.path.exists(os.path.join(self.output_directory, 'factors{}.dat'.format(zone))):
            #    os.remove(os.path.join(self.output_directory, 'factors{}.dat'.format(zone)))
                
            self.run_ppk2fac(self.ppk2fac_exe, 
                             instruct_fname='ppk2fac{}.in'.format(zone))
            self.write_fac2real_instruct('factors{}.dat'.format(zone),
                                    'points{}.pts'.format(zone), 
                                    1E-6, 1E6, 'values{}.ref'.format(zone),
                                    instruct_fname='fac2real{}.in'.format(zone))

    def run_pyfac2real_by_zones(self, zones):
        '''
        Function to run fac2real via the python implementation and based on 
        number of zones.
        
        This should only be used after running setup_pilot_points_by_zones
        '''      
        for zone in range(zones):
            count = 0
            while not os.path.exists(os.path.join(self.output_directory, "factors{}.dat".format(zone))):
                print('Factors file not found ...waiting 1 more second')
                count += 1
                time.sleep(1)
                if count == 5:
                    break
            self.run_pyfac2real(pp_fname="points{}.pts".format(zone),
                              factors_fname="factors{}.dat".format(zone), 
                              out_fname="values{}.ref".format(zone),
                              upper_lim=1.0e+6, lower_lim=1.0e-6)        

            
        self._create_mesh3D_array_from_values(self.zone_array, values_dir=self.output_directory, 
                                              out_dir=self.output_directory)

    def update_pilot_points_files_by_zones(self, zones, new_values_dict):
        
        for zone in range(zones):
            self.update_pilot_points_file_by_zone(new_values_dict[zone], zone,
                                points_fname='points{}.pts'.format(zone))
        
# End class PilotPoints
        
def mesh3DToVtk(mesh_array, grid_info, val_array, val_name, out_path, vtk_out):
    '''
    Function to write the mesh array 
    '''
    from HydroModelBuilder.GISInterface.GDALInterface import array2Vtk

    mesh, zone_matrix = mesh_array[0], mesh_array[1]
    array2Vtk.build_vtk_from_array(grid_info, np.fliplr(mesh), ["z_elev"], 
                                   [np.fliplr(mesh)], ["zone", val_name], 
                                   [np.fliplr(zone_matrix), np.fliplr(val_array)], 
                                   out_path, vtk_out)
    
if __name__ == "__main__":

    from HydroModelBuilder.GWModelManager import GWModelManager

    resolution = 5000
    zone_map = {1: 'qa', 2: 'utb', 3: 'utqa', 4: 'utam', 5: 'utaf', 6: 'lta', 
                7: 'bse'}
    HGU_map = {'bse':'Bedrock', 'utb':'Newer Volcanics Basalts', 
               'utaf':'Calivil', 'lta':'Renmark', 
               'qa':'Coonambidgal Sands', 'utqa':'Shepparton Sands',
               'utam':'Loxton-Parilla Sands'}

    
    # Get an example mesh from a previous model build
    MM = GWModelManager()
    model_folder = r'C:/Workspace/part0075/MDB modelling/testbox/00_Campaspe_Cascade/01_steady_state/structured_model_grid_{}m'.format(resolution)
    MM.load_GW_model(os.path.join(model_folder, r"01_steady_state_packaged.pkl"))
    name = MM.GW_build.keys()[0]
    mesh_array = MM.GW_build[name].model_mesh3D
    zones = len(np.unique(mesh_array[1])) - 1
    cell_centers = MM.GW_build[name].model_mesh_centroids
    model_boundary = MM.GW_build[name].model_boundary

    # Generate pilot points from existing mesh with zones and cell centers
    # allow for assignment of pilot points at every nth active zone cell using
    # skip function.
    
    pp = PilotPoints(output_directory=model_folder)

    a = pp.generate_points_from_mesh(mesh_array, cell_centers, 
                   skip=[0,  0, 3, 0, 2, 3, 3], 
            skip_active=[3, 20, 0, 4, 0, 0, 0],
            zone_prop_dict={0:30.0, 1:2.0, 2:2.0, 3:50.0, 4:45.0, 5:30.0, 6:10.0},
            add_noise=True
            )

    pp.write_settings_fig()
    pp.write_grid_spec(mesh_array, model_boundary, delc=resolution, delr=resolution)
    pp.write_struct_file(mesh_array, nugget=0.0, 
                      transform='log',numvariogram=1, variogram=0.15, 
                      vartype=2, bearing=0.0, a=20000.0, anisotropy=1.0)

    search_radius = [30000, 20000, 20000, 20000, 20000, 20000, 20000]

    pp.setup_pilot_points_by_zones(mesh_array, zones, search_radius)    

    pp.run_pyfac2real_by_zones(zones)
#    for zone in range(zones):
#        print('There are {0} points in zone: {1}'.format(pp.num_ppoints_by_zone[zone], zone))
#        
##        run_fac2real(fac2real_exe, instruct_fname='fac2real{}.in'.format(zone))
#        pp.run_pyfac2real(pp_fname="points{}.pts".format(zone),
#                          factors_fname="factors{}.dat".format(zone), 
#                          out_fname="values{}.ref".format(zone),
#                          upper_lim=1.0e+6, lower_lim=1.0e-6)


    hk = pp.val_array

    nrow, ncol = mesh_array[1][0].shape[0], mesh_array[1][0].shape[1]    
    delc, delr = [resolution]*2
    x0, y0 = model_boundary[0], model_boundary[3]
    grid_info = [ncol, nrow, delc, delr, x0, y0]
    mesh3DToVtk(mesh_array, grid_info, hk, 'hk', '', 'hk')
        
    ###########################################################################
    from HydroModelBuilder.ModelInterface.flopyInterface import flopyInterface
    data_folder = r"C:/Workspace/part0075/MDB modelling/testbox/PEST5000/master"
    modflow_model = flopyInterface.ModflowModel(MM.GW_build[name], data_folder=data_folder)
    modflow_model.buildMODFLOW()
    import flopy
    import matplotlib.pyplot as plt
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
        array = modelmap.plot_array(vals, masked_values=[1E6], alpha=1.0)#, 
                                    #vmin=vmin, vmax=vmax)
        start, end = ax.get_xlim()
        start = start // 1000 * 1000 + 1000
        end = end // 1000 * 1000 - 1000
        ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
        plt.scatter([x[0] for x in b[key]], [x[1] for x in b[key]], facecolors='none', color='black', alpha=0.5)
        cbar_ax = fig.add_axes() #[0.67, 0.055, 0.01, 0.42])
        fig.colorbar(array, cax=cbar_ax)        
    
    fig.subplots_adjust(left=0.01, right=0.95, bottom=0.05, top=0.95, wspace=0.1, hspace=0.12)