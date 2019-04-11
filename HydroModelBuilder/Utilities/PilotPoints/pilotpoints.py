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


class PilotPointsLinear(object):
    """This class is used to handle linear "pilot points" in a simple fashion"""

    def __init__(self):
        pass

    def set_uniform_points(self, length, num_points):
        '''
        Generate a list of length num_points, that has even increments length / num_points
        param: length, float of length
        param: num_points, integer

        '''
        self.length = length
        self.num_points = num_points
        self.points = [length / float(num_points - 1) * x for x in range(num_points)]

    def interpolate_unknown_points_from_df_col(self, df_col, points_vals):
        """TODO: Docs

        :param df_col:
        :param points_vals:
        """
        return np.interp(df_col.tolist(), self.points, points_vals)


class PilotPoints(object):
    """TODO: Docs"""

    def __init__(self, output_directory=None, additional_name=""):
        self.min_allowable_points_separation = 0.0
        self.grid_spec_fname = 'grid.spc'
        self.points_fname = 'points.pts'
        self.zone_fname = 'zone.inf'
        self.struct_fname = 'struct.dat'
        self.ppk2fac_exe = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        'ppk2fac.exe')
        #self.ppk2fac_exe = os.path.join(os.getcwd(), 'ppk2fac.exe')
        self.fac2real_exe = 'fac2real.exe'
        self.ppcov_exe = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      'ppcov.exe')
        self.output_directory = os.path.join(output_directory, '{}_pilot_points'.format(additional_name))
        if not os.path.exists(self.output_directory):
            os.mkdir(self.output_directory)

        self.num_ppoints_by_zone = {}

    def write_settings_fig(self, colrow='yes', date=r'mm/dd/yyyy'):
        """Write the "settings.fig" file that is required to be present when running ppk2fac

        :param colrow:  (Default value = 'yes')
        :param date:  (Default value = r'mm/dd/yyyy')
        """
        out_fname = os.path.join(self.output_directory, 'settings.fig')
        with open(out_fname, 'w') as f:
            f.write('colrow={} \n'.format(colrow))
            f.write('date={}'.format(date))

    @staticmethod
    def _zone_array2layers(zone_array):
        """Generate masked layer arrays for each zone.

        :param zone_array:
        """
        #zones = [x + 1 for x in range(len(np.unique(zone_array)) - 1)]
        zones = np.unique(zone_array).astype(int).tolist()
        zones.remove(-1)
        layers = zone_array.shape[0]
        zone_mask2D = {}
        for index, zone in enumerate(zones):
            zone_mask = np.ma.masked_array(zone_array,
                                           zone_array == zone).mask

            zone_mask2D[zone] = np.full_like(zone_mask[0], False, dtype=bool)
            for layer in range(layers):
                zone_mask2D[zone] |= zone_mask[layer]

        return zone_mask2D

    def generate_points_from_mesh(self, mesh_array, cell_centers, zone_prop_dict=None,
                                  skip=0, skip_active=0, add_noise=False, zones=None,
                                  force_single_point=False):
        """Generate pilot points based in the zonal array in the mesh array object.

        :param mesh_array:
        :param cell_centers:
        :param zone_prop_dict:  (Default value = None)
        :param skip:  (Default value = 0)
        :param skip_active:  (Default value = 0)
        :param add_noise:  (Default value = False)

        :returns: array, points zone array and initial values.
        """
        # First get the number of zones and layers from the mesh_array object
        if zones is None:
            zones = [x + 1 for x in range(len(np.unique(mesh_array[1])) - 1)]
        
        points_dict = {}
        points_zone_dict = {}
        points_val_dict = {}

        zone_mask2D = self._zone_array2layers(mesh_array[1])
        zone_active = []
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

            points_dict[zone] = []
            points_zone_dict[zone] = []
            points_val_dict[zone] = []

            row_points = 0
            col_points = 0
            points_active = 0
            for row in range(zone_mask2D[zone].shape[0]):
                if skip != 0:
                    row_points += 1
                    if row_points > 1 and row_points < skip_n:
                        continue
                    elif row_points == skip_n:
                        row_points = 0
                        continue
                for col in range(zone_mask2D[zone].shape[1]):
                    if skip_n != 0:
                        col_points += 1
                        if col_points > 1 and col_points < skip_n:
                            continue
                        elif col_points == skip_n:
                            col_points = 0
                            continue

                    if zone_mask2D[zone][row][col]:
                        if skip_active_n != 0:
                            points_active += 1
                            if points_active > 1 and points_active < skip_active_n:
                                continue
                            elif points_active == skip_active_n:
                                points_active = 0
                                continue

                        points_dict[zone] += [(cell_centers[0][row][col],
                                                cell_centers[1][row][col])]

                        points_zone_dict[zone] += [zone]
                        if zone_prop_dict is not None:
                            x = zone_prop_dict[zone]
                            if add_noise:
                                points_val_dict[zone] += [x + 0.1 * x * np.random.standard_normal()]
                            else:
                                points_val_dict[zone] += [x]
                        else:
                            if add_noise:
                                points_val_dict[zone] += [1.0 + np.random.lognormal(-0.5, 0.2)]
                            else:
                                points_val_dict[zone] += [1.0]
            
            # Test that points were derived and if not pop it from the list but with warning
            if points_dict[zone] == []:
                points_dict.pop(zone)
                points_zone_dict.pop(zone)
                points_val_dict.pop(zone)
                print("WARNING! No points were determined for zone: {}, removing".format(zone))
            else:
                zone_active += [zone]

        self.points_dict = points_dict
        self.points_zone_dict = points_zone_dict
        self.points_val_dict = points_val_dict
        self.zone_mask2D = zone_mask2D
        return points_dict, points_zone_dict, points_val_dict, zone_mask2D, zone_active

    def write_grid_spec(self, mesh_array, model_boundary,
                        grid_spec_fname='grid.spc', delc=None, delr=None):
        """Write the grid specification file using a mesh array for use in the PPK2FAC utility.

        This only works for model grids that are evenly spaced along the rows and columns

        ** Future implementation should allow for any type of array which will need
           to be inherited from the model_mesh object in GWModelBuilder

        :param mesh_array:
        :param model_boundary:
        :param grid_spec_fname:  (Default value = 'grid.spc')
        :param delc:  (Default value = None)
        :param delr:  (Default value = None)
        """
        layer = 0
        row, col = mesh_array[1][layer].shape[0], mesh_array[1][layer].shape[1]
        with open(os.path.join(self.output_directory, grid_spec_fname), 'w') as f:
            f.write('{0} {1}\n'.format(row, col))
            f.write('{0} {1} {2} \n'.format(model_boundary[0], model_boundary[3], 0.0))
            f.write(" ".join(['{}'.format(delc) for x in range(col)]) + '\n')
            f.write(" ".join(['{}'.format(delr) for x in range(row)]) + '\n')

    def write_pilot_points_file(self, pilot_points, pp_zones, initial_values,
                                points_fname='points.pts', prefix=None):
        """Write the pilot point files containing name, easting, northing, zone and value for each pilot point.

        :param pilot_points:
        :param pp_zones:
        :param initial_values:
        :param points_fname:  (Default value = 'points.pts')
        :param prefix:  (Default value = None)
        """
        with open(os.path.join(self.output_directory, points_fname), 'w') as f:
            for index, point in enumerate(pilot_points):
                if prefix is not None:
                    f.write('{0}{1} {2} {3} {4} {5} \n'.format(prefix, index, point[0],
                                                               point[1], pp_zones[index] + 1, initial_values[index]))
                else:
                    f.write('pp{0} {1} {2} {3} {4} \n'.format(index, point[0],
                                                              point[1], pp_zones[index] + 1, initial_values[index]))

    def update_pilot_points_file_by_zone(self, new_values, zone,
                                         points_fname='points.pts', prefix=None):
        """Write the pilot point files containing name, easting, northing, zone and value for each pilot point

        :param new_values:
        :param zone:
        :param points_fname:  (Default value = 'points.pts')
        :param prefix:  (Default value = None)
        """
        with open(os.path.join(self.output_directory, points_fname), 'r') as f:
            prefix = f.readlines()[0].split()[0][:-1]

        with open(os.path.join(self.output_directory, points_fname), 'w') as f:
            for index, point in enumerate(self.points_dict[zone]):
                if prefix is not None:
                    f.write('{0}{1} {2} {3} {4} {5} \n'.format(prefix, index, point[0],
                                                               point[1], self.points_zone_dict[zone][index] + 1, new_values[index]))
                else:
                    f.write('pp{0} {1} {2} {3} {4} \n'.format(index, point[0],
                                                              point[1], self.points_zone_dict[zone][index] + 1, new_values[index]))

    def write_zone_file(self, zone_array, zone_fname='zone.inf'):
        """Write the zone array file for use in ppk2fac based on an integer array of cell zones.

        :param zone_array:
        :param zone_fname:  (Default value = 'zone.inf')
        """
        layer = zone_array.keys()[0]
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
                          transform='log', numvariogram=1, variogram=0.15,
                          vartype=2, bearing=0.0, a=500.0, anisotropy=1.0):
        """Write the structure file for use in ppk2fac using default values as specified
        in the arguments to the function, which can then be manually modified through "struct_file_out"

        :param mesh_array:
        :param struct_fname:  (Default value = 'struct.dat')
        :param nugget:  (Default value = 0.0)
        :param transform:  (Default value = 'log')
        :param numvariogram:  (Default value = 1)
        :param variogram:  (Default value = 0.15)
        :param vartype:  (Default value = 2)
        :param bearing:  (Default value = 0.0)
        :param a:  (Default value = 500.0)
        :param anisotropy:  (Default value = 1.0)
        """
        zones = np.unique(mesh_array[1]).astype(int).tolist()
        zones.remove(-1)
        with open(os.path.join(self.output_directory, struct_fname), 'w') as f:
            for zone in zones:
                f.write('STRUCTURE structure{}\n'.format(zone))
                f.write('  NUGGET {}\n'.format(nugget))
                f.write('  TRANSFORM {}\n'.format(transform))
                f.write('  NUMVARIOGRAM {}\n'.format(numvariogram))
                f.write('  VARIOGRAM vario{0} {1}\n'.format(zone, variogram))
                f.write('END STRUCTURE\n')
                f.write('\n')
            for zone in zones:
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
        """Write the instructions to be used with executing ppk2fac to allow batch execution.

        :param mesh_array:
        :param grid_spec_fname:
        :param points_fname:
        :param struct_fname:
        :param zone_fname:  (Default value = None)
        :param instruct_fname:  (Default value = 'ppk2fac.in')
        :param min_allowable_points_separation:  (Default value = 0.0)
        :param kriging:  (Default value = 'o')
        :param search_radius:  (Default value = 10000)
        :param min_pilot_points:  (Default value = 1)
        :param max_pilot_points:  (Default value = 5)
        :param factor_fname:  (Default value = 'factors.dat')
        :param formatted:  (Default value = 'f')
        :param sd_fname:  (Default value = 'sd.ref')
        :param reg_fname:  (Default value = 'reg.dat')
        """
        zones = np.unique(mesh_array[1]).astype(int).tolist()
        zones.remove(-1)
        with open(os.path.join(self.output_directory, instruct_fname), 'w') as f:
            f.write('{} \n'.format(grid_spec_fname))
            f.write('{} \n'.format(points_fname))
            f.write('{} \n'.format(min_allowable_points_separation))
            if zone_fname is not None:
                f.write('{} \n'.format(zone_fname))
            else:
                f.write(' \n')
            # end if
            f.write('{} \n'.format(struct_fname))
            f.write('\n')  # First zone represents inactive cells so skipped
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
        """Run ppk2fac using generated input commands from write_ppk2fac_instruct

        :param ppk2fac_exe:
        :param instruct_fname:  (Default value = None)
        """
        # Change into location where input files for ppk2fac are located
        cwd = os.getcwd()
        os.chdir(self.output_directory)

        # Run ppk2fac
        if instruct_fname is not None:
            subprocess.Popen([ppk2fac_exe, '<', instruct_fname], shell=True)
        else:
            command = ppk2fac_exe
            subprocess.call(command)
        # end if

        # Return to current working directory
        os.chdir(cwd)

    def write_fac2real_instruct(self, factors_file, points_file, lower_lim,
                                upper_lim, out_file,
                                instruct_fname='fac2real.in'):
        """Write the instructions to be used with executing fac2real
        - Assumes formatted file for factors file
        - Assumes use of singular value rather than array
        - Sets value for elements to which no interpolation takes place

        :param factors_file:
        :param points_file:
        :param lower_lim:
        :param upper_lim:
        :param out_file:
        :param instruct_fname:  (Default value = 'fac2real.in')
        """
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
        """ run fac2real using generated input commands from write_fac2real_instruct

        :param fac2real_exe:
        :param instruct_fname:  (Default value = None)
        """
        # Change into location where input files for ppk2fac are located
        cwd = os.getcwd()
        os.chdir(self.output_directory)

        # Run fac2real
        if instruct_fname is not None:
            subprocess.call([fac2real_exe, '<', instruct_fname], shell=True)
        else:
            command = fac2real_exe
            subprocess.call(command)
        # end if

        # Return to current working directory
        os.chdir(cwd)

    def run_pyfac2real(self, pp_fname='points.dat', factors_fname='factors.dat',
                       out_fname='values.ref', lower_lim=1.0E-6, upper_lim=1.0E+6):
        """Run Python implementation of fac2real using the implementation in pyEmu by Jeremy White, see:
        https://github.com/jtwhite79/pyemu/blob/master/pyemu/utils/gw_utils.py

        :param pp_fname:  (Default value = 'points.dat')
        :param factors_fname:  (Default value = 'factors.dat')
        :param out_fname:  (Default value = 'values.ref')
        :param lower_lim:  (Default value = 1.0E-6)
        :param upper_lim:  (Default value = 1.0E+6)
        """
        import pyemu

        pyemu.utils.fac2real(str(os.path.join(self.output_directory, pp_fname)),
                             str(os.path.join(self.output_directory, factors_fname)),
                             out_file=os.path.join(self.output_directory, out_fname),
                             upper_lim=upper_lim, lower_lim=lower_lim)

    def write_ppcov_instruct(self, pp_fname, struct_fname, structure_name,
                             min_allowable_points_separation=0.0,
                             pp_param_prefix=None, out_fname='points.mat',
                             instruct_fname='ppcov.in'):
        """Write the instructions to be used with executing "ppcov" for generating
        covariance matrix for pilot points.

        - Assumes formatted file for factors file
        - Assumes use of singular value rather than array
        - Sets value for elements to which no interpolation takes place

        :param pp_fname:
        :param struct_fname:
        :param structure_name:
        :param min_allowable_points_separation:  (Default value = 0.0)
        :param pp_param_prefix:  (Default value = None)
        :param out_fname:  (Default value = 'points.mat')
        :param instruct_fname:  (Default value = 'ppcov.in')
        """

        with open(os.path.join(self.output_directory, instruct_fname), 'w') as f:
            f.write('{} \n'.format(pp_fname))
            f.write('{} \n'.format(min_allowable_points_separation))
            f.write('{} \n'.format(struct_fname))
            f.write('{} \n'.format(structure_name))
            f.write('{} \n'.format(out_fname))
            if pp_param_prefix is not None:
                f.write('{} \n'.format(pp_param_prefix))
            else:
                f.write('\n')

    def run_ppcov(self, ppcov_exe, instruct_fname=None):
        """Run ppcov using generated input commands from write_ppcov_instruct

        :param ppcov_exe:
        :param instruct_fname:  (Default value = None)
        """

        # Change into location where input files for ppk2fac are located
        cwd = os.getcwd()
        os.chdir(self.output_directory)
        # Run ppcov
        if instruct_fname is not None:
            subprocess.call([ppcov_exe, '<', instruct_fname], shell=True)
        else:
            # Note that this will bring up an interactive window requiring
            # inputs
            command = ppcov_exe
            subprocess.call(command)
        # end if

        # Return to current working directory
        os.chdir(cwd)

    def _create_mesh3D_array_from_values(self, zone_array, values_dir='', out_dir=''):
        """Read in all of the different values arrays for each layer and
        construct the 3D array pertaining to pilot point values as defined by their zones

        :param zone_array:
        :param values_dir:  (Default value = '')
        :param out_dir:  (Default value = '')
        """
        val_array = np.zeros_like(zone_array)
        #zones = [x + 1 for x in range(len(np.unique(zone_array)) - 1)]
        zones = self.points_dict.keys()         
        for index, zone in enumerate(zones):
            fname = os.path.join(self.output_directory, 'values{}.ref'.format(zone))
            with open(fname, 'r') as f:
                lines = f.readlines()
                vals = []
                for ind, line in enumerate(lines):
                    vals += [float(x) for x in line.strip('\n').split()]

            vals = np.array(vals)
            vals = np.reshape(vals, zone_array[0].shape)
            for layer in range(zone_array.shape[0]):
                val_array[layer][zone_array[layer] == [zone]] = vals[zone_array[layer] == [zone]]

        self.val_array = val_array

        # return val_array

    def save_mesh3D_array(self, filename='val_array', data_format='binary'):
        """Save the created 3D mesh array to file with filename as 'filename' 
        and as either 'data_format' as 'binary' or 'ascii'

        :param filename:  (Default value = 'val_array')
        :param data_format:  (Default value = 'binary')
        """
        if data_format == 'ascii':
            np.savetxt(filename, self.val_array)
        elif data_format == 'binary':
            np.save(filename, self.val_array)
        else:
            print 'Data format not recognised, use "binary" or "ascii"'
        # end if

    def setup_pilot_points_by_zones(self, mesh_array, zones, search_radius,
                                    verbose=True, prefixes=None):
        """Set up pilot points files, zone files, ppk2fac instructions, and to run ppk2fac
        to create the factors files required for fac2real.

        This is based on creating a pilot point grid per zone in the mesh.

        :param mesh_array:
        :param zones:
        :param search_radius:
        :param verbose:  (Default value = True)
        :param prefixes:  (Default value = None)
        """
        if type(zones) == int:
            zones = range(zones)
        # end if
        
        for index, zone in enumerate(zones):
            if prefixes is not None:
                prefix = prefixes[index]
            else:
                prefix = None
            print prefix
            self.write_pilot_points_file(self.points_dict[zone],
                                         self.points_zone_dict[zone],
                                         self.points_val_dict[zone],
                                         points_fname='points{}.pts'.format(zone),
                                         prefix=prefix)

            self.num_ppoints_by_zone[zone] = len(self.points_dict[zone])

            self.write_zone_file(self.zone_mask2D,
                                 zone_fname='zone{}.inf'.format(zone))

            self.write_ppk2fac_instruct(mesh_array, 'grid.spc',
                                        'points{}.pts'.format(zone),
                                        zone_fname='zone{}.inf'.format(zone),
                                        instruct_fname='ppk2fac{}.in'.format(zone),
                                        struct_fname='struct.dat',
                                        factor_fname='factors{}.dat'.format(zone),
                                        sd_fname='sd{}.ref'.format(zone),
                                        reg_fname='reg{}.dat'.format(zone),
                                        min_allowable_points_separation=0.0, kriging='o',
                                        search_radius=search_radius[index],
                                        min_pilot_points=1,
                                        max_pilot_points=50)

            # if os.path.exists(os.path.join(self.output_directory, 'factors{}.dat'.format(zone))):
            #    os.remove(os.path.join(self.output_directory, 'factors{}.dat'.format(zone)))

        for zone in zones:
            self.run_ppk2fac(self.ppk2fac_exe,
                             instruct_fname='ppk2fac{}.in'.format(zone))

            self.write_fac2real_instruct('factors{}.dat'.format(zone),
                                         'points{}.pts'.format(zone),
                                         1E-6, 1E6, 'values{}.ref'.format(zone),
                                         instruct_fname='fac2real{}.in'.format(zone))

    def generate_cov_mat_by_zones(self, zones):
        """Set up ppcov instructions, then run ppcov to create the covmat files
        required for building the covariance matrix.

        This is based on creating a pilot point grid per zone in the mesh.

        :param zones: int of number of zones or list of zone ids
        """
        if type(zones) == int:
            zones = range(zones)
        # end if

        for zone in zones:
            self.write_ppcov_instruct('points{}.pts'.format(zone),
                                      'struct.dat',
                                      'structure{}'.format(zone),
                                      out_fname='points{}.mat'.format(zone),
                                      instruct_fname='ppcov{}.in'.format(zone)
                                      )

        for zone in zones:
            self.run_ppcov(self.ppcov_exe, instruct_fname='ppcov{}.in'.format(zone))

    def run_pyfac2real_by_zones(self, zones):
        """Run fac2real via the python implementation and based on number of zones.

        This should only be used after running setup_pilot_points_by_zones

        :param zones: int of number of zones or list of zone ids
        """
        
        if type(zones) == int:
            zones = range(1, zones+1)
        # end if
        
        for zone in zones:
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
        """TODO: Docs

        :param zones:
        :param new_values_dict:
        """
        for zone in range(zones):
            self.update_pilot_points_file_by_zone(new_values_dict[zone], zone+1,
                                                  points_fname='points{}.pts'.format(zone+1))

# End PilotPoints()


def mesh3DToVtk(mesh_array, grid_info, val_array, val_name, out_path, vtk_out):
    """Write the mesh array.

    :param mesh_array:
    :param grid_info:
    :param val_array:
    :param val_name:
    :param out_path:
    :param vtk_out:
    """
    from HydroModelBuilder.GISInterface.GDALInterface import array2Vtk

    mesh, zone_matrix = mesh_array[0], mesh_array[1]
    array2Vtk.build_vtk_from_array(grid_info, np.fliplr(mesh), ["z_elev"],
                                   [np.fliplr(mesh)], ["zone", val_name],
                                   [np.fliplr(zone_matrix), np.fliplr(val_array)],
                                   out_path, vtk_out)


if __name__ == "__main__":
    # Shifted testing script to ../../tests folder
    pass
