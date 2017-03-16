"""
Interface for PEST software to generate necessary files
Based on code of Etienne Bresciani (applied for mflab in Matlab) translated to 
python and modified as well to work with the GWModelBuilder class.
"""

import os
import sys
import datetime
import csv
import numpy as np
import pandas as pd

class PESTInterface(object):

    def __init__(self, name=None, directory=None, csv_copy=False, 
                 excel_copy=False, params=None, obs=None, obs_grp=None, 
                 models_ID=None, predictive_analysis=False):
        self.PEST_data = {}
        if name == None:
            self.name = 'default' 
        else:
            self.name = name
        #end if
        if directory == None:
            self.directory = ''
        else:
            self.directory = directory
        #end if
            
        self.csv_copy = csv_copy
        self.excel_copy = excel_copy
        if params:
            self.params = params
        else:
            self.params = {}
        #end if
        if obs:
            self.obs = obs
        else:
            self.obs = {}
        #end if

        if obs_grp:
            self.obs_grp = obs_grp
        else:
            self.obs_grp = {}
        #end if
            
        #self.obs = obs

        self.predictive_analysis = predictive_analysis
        
        self.models_ID = models_ID        
        
        self.PEST_data['PESTcon'] = {}
        self.PESTcon()
        self.PEST_data['PESTpar'] = self.PESTpar(params=self.params)
        
        #self.PEST_data['PESTpgp'] = {}
        #self.PESTpgp()
        self.PEST_data['PESTobs'] = self.PESTobs(obs=self.obs, obs_grp=self.obs_grp)
        
    def PESTcon(self):
        """
        Control data for the *.pst file of pest, which takes the form:

        * control data
        RSTFLE PESTMODE
        NPAR NOBS NPARGP NPRIOR NOBSGP [MAXCOMPDIM] [DERZEROLIM]
        NTPLFLE NINSFLE PRECIS DPOINT [NUMCOM JACFILE MESSFILE] [OBSREREF]
        RLAMBDA1 RLAMFAC PHIRATSUF PHIREDLAM NUMLAM [JACUPDATE] [LAMFORGIVE] [DERFORGIVE]
        RELPARMAX FACPARMAX FACORIG [IBOUNDSTICK UPVECBEND] [ABSPARMAX]
        PHIREDSWH [NOPTSWITCH] [SPLITSWH] [DOAUI] [DOSENREUSE] [BOUNDSCALE]
        NOPTMAX PHIREDSTP NPHISTP NPHINORED RELPARSTP NRELPAR [PHISTOPTHRESH] [LASTRUN] [PHIABANDON]
        ICOV ICOR IEIG [IRES] [JCOSAVE] [VERBOSEREC] [JCOSAVEITN] [REISAVEITN] [PARSAVEITN] [PARSAVERUN]

        Where (brief description of options, but details can be accessed from
        the PEST user manual):
        
        *** 2nd line ***
        RSTFLE = "restart" or "norestart", with the former allowing PEST to 
        recommence if the PEST run is halted [string]
        
        PESTMODE = "estimation", "prediction", "regularisation" or "pareto" [string]

        *** 3rd line ***
        NPAR  = total number of parameters [integer] (self generated)
        
        NOBS = the total number of observations [integer] (self generated)
        
        NPARGP = number of parameter groups [integer] (self generated)
        
        NPRIOR = number of articles of prior information included in the parameter
        estimation process [integer] (self generated)
        
        NOBSGP = number of observation groups [integer] (self generated)

        MAXCOMPDIM (optional) = activates compressed internal storage of the Jacobian matrix [integer]

        DERZEROLIM (optional) = threshold for considering element of the Jacobian matrix zero [float]         
        
        *** 4th line ***
        NTPLFLE = number of template files [integer]
       
        NINSFLE = number of instruction files [integer]
        
        PRECIS = "single" or "double", precision with which PEST writes parameters to input file [string]

        DPOINT = "point" or "nopoint", allows ignoring of decimal points if the latter option is chosen [string]         
        
        NUMCOM, JACFILE and MESSFILE (optional) = the manner in which PEST can obtain derivatives 
          directly from the model, typically set as 1,0,0 [3 * integer]       
        
        OBSREREF (optional) = "obsreref" or "noobsreref" for observation re-referencing [string]
        
        *** 5th line ***
        RLAMBDA1 = initial Marquardt lambda [real]

        RLAMFAC = factor for adjusting the Marquardt lambda, set as >1.0 or <-1.0 [float]                
        
        PHIRATSUF = stands for "phi ratio sufficient", real variable, 0.3 is mostly appropriate,  [float]
        
        ...
        
        *** end file ***
        
        NOTE: Run pestchek to insure that all control data has been entered appropriately
        """
        
        control_data = {'RSTFLE': 'restart',
                        'PESTMODE': 'estimation',
                        'PRECIS': 'single',
                        'DPOINT': 'point',
                        'RLAMBDA1': 20,
                        'RLAMFAC': -3,
                        'PHIRATSUF': 0.3,
                        'PHIREDLAM': 1.00E-02,
                        'NUMLAM': 7,
                        'JACUPDATE': 999,
                        'LAMFORGIVE': 'lamforgive',
                        'DERFORGIVE': 'derforgive',
                        'RELPARMAX': 10,
                        'FACPARMAX': 10,
                        'FACORIG': 1.00E-03,
                        'PHIREDSHW': 0.1,
                        'NOPTSWITCH': 1,
                        'BOUNDSCALE': 'noboundscale', #'boundscale',
                        'NOPTMAX':	25,
                        'PHIREDSTP':	0.01,
                        'NPHISTP':	5,
                        'NPHINORED': 5,
                        'RELPARSTP': 0.01,
                        'NRELPAR': 3,
                        'ICOV': 0,
                        'ICOR': 0,
                        'IEIG': 0,
                        'IRES': 0,
                        'JCOSAVE': 'jcosave',
                        'VERBOSEREC': 'verboserec',
                        'JCOSAVEITN': 'jcosaveitn',
                        'REISAVEITN': 'reisaveitn',
                        'PARSAVEITN': 'parsaveitn'
                        }        
                        
        singular_value_decomposition = {'SVDMODE': 1,
                                        'MAXSING': 0, # Number of parameters
                                        'EIGTHRESH': 5.00E-07,
                                        'EIGWRITE': 1,
                                        }


        predictive_analysis = {'NPREDMAXMIN':	1,
                               'PD0': 0,
                               'PD1': 0.00E+00,
                               'PD2': 0,
                               'ABSPREDLAM': 0,
                               'RELPREDLAM': 0.01,
                               'INITSCHFAC': 1,
                               'MULSCHFAC': 2,
                               'NSEARCH': 6,
                               'ABSPREDSWH': 0,
                               'RELPREDSWH': 0.05,
                               'NPREDNORED': 4,
                               'ABSPREDSTP': 0,
                               'RELPREDSTP': 0.005,
                               'NPREDSTP': 4
                               }

        self.PEST_data['PESTcon'] = {}
        self.PEST_data['PESTcon']['control_data'] = control_data
        self.PEST_data['PESTcon']['singular_value_decomposition'] = singular_value_decomposition
        self.PEST_data['PESTcon']['predictive_analysis'] = predictive_analysis

    def PESTpar(self, params=None):
        """ 
        Generate the PEST parameter file
        
        For column 'PARNAME' the maximum length is 12 characters
        PARTRANS options = ['log', 'fixed', 'tied']
        PARCHGLIM options = ['factor', 'relative']
        
        """
        header = ['PARNAME', 'PARTRANS', 'PARCHGLIM', 'PARVAL1', 'PARLBND', 'PARUBND', 'PARGP', 'SCALE', 'OFFSET', 'PARTIED', 'models', 'unit', 'comment']
        num_param = len(params.keys())
        PESTpar = pd.DataFrame(columns=header, index=self.params.keys())
        PESTpar['PARNAME'] = params.keys()
        PESTpar['PARTRANS'] = [x['PARTRANS'] if 'PARTRANS' in x.keys() else 'log' for x in params.values()] #['log'] * num_param 
        PESTpar['PARCHGLIM'] = [x['PARCHGLIM'] if 'PARCHGLIM' in x.keys() else 'factor' for x in params.values()] # ['factor'] * num_param         
        PESTpar['PARVAL1'] = [x['PARVAL1'] for x in params.values()]
        PESTpar['PARLBND'] = [x['PARLBND'] if 'PARLBND' in x.keys() else x['PARVAL1'] * 0.9 for x in params.values()] #[0] * num_param        
        PESTpar['PARUBND'] = [x['PARUBND'] if 'PARUBND' in x.keys() else x['PARVAL1'] * 1.1 for x in params.values()] #[0] * num_param
        PESTpar['PARGP'] = [x['PARGP'] if 'PARGP' in x.keys() else 'default' for x in params.values()] #['default'] * num_param
        PESTpar['SCALE'] = [x['SCALE'] if 'SCALE' in x.keys() else 1.0 for x in params.values()] #[1.0] * num_param
        PESTpar['OFFSET'] = [x['OFFSET'] if 'OFFSET' in x.keys() else 0.0 for x in params.values()] #[0.0] * num_param
        PESTpar['PARTIED'] = [''] * num_param
        PESTpar['models'] = [self.models_ID[0]] * num_param
        PESTpar['unit'] = ['-'] * num_param
        
        if self.csv_copy:
            self._createCSVcopyFromDataFrame(PESTpar, 'PESTpar')
        #end if
        return PESTpar

    def genPESTpgp(self):
        header = ['PARGPNME', 'INCTYP', 'DERINC', 'DERINCLB', 'FORCEN', 'DERINCMUL', 'DERMTHD']
        INCTYPdefault = 'relative'
        DERINCdefault = 0.01
        DERINCLBdefault = 0.0001
        FORCENdefault = 'switch'
        DERINCMULdefault = 1.5
        DERMTHDdefault = 'parabolic'        
        # Get number of parameter groups from PESTpar        
        parameter_groups = self.PEST_data['PESTpar']['PARGP'].unique()        
        num_parameter_groups = len(parameter_groups)

        PESTpgp = pd.DataFrame(columns=header, index=parameter_groups)
        PESTpgp['PARGPNME'] = parameter_groups
        PESTpgp['INCTYP'] = [INCTYPdefault] * num_parameter_groups
        PESTpgp['DERINC'] =  [DERINCdefault] * num_parameter_groups
        PESTpgp['DERINCLB'] =  [DERINCLBdefault] * num_parameter_groups
        PESTpgp['FORCEN'] =  [FORCENdefault] * num_parameter_groups
        PESTpgp['DERINCMUL'] =  [DERINCMULdefault] * num_parameter_groups
        PESTpgp['DERMTHD'] =  [DERMTHDdefault] * num_parameter_groups       

        if self.csv_copy:
            self._createCSVcopyFromDataFrame(PESTpgp, 'PESTpgp')
        #end if
        self.PEST_data['PESTpgp'] = PESTpgp
        
    def PESTobs(self, obs=None, obs_grp=None):
        
        for index, key in enumerate(obs_grp.keys()):
            obs_grp_ts = obs_grp[key]['time_series'].copy()
            # Filter out null observations that were out of date range or 
            # didn't map to the model mesh
            obs_grp_ts = obs_grp_ts[obs_grp_ts['obs_map'] != 'null']
            if obs_grp[key]['by_zone']:
                obs_grp_ts['OBGNME'] = obs_grp_ts['zone']
            else:
                obs_grp_ts['OBGNME'] = key 
            # end if            
            obs_grp_ts['WEIGHT'] = obs_grp[key]['weights']
            obs_grp_ts['model'] = self.models_ID[0]
            obs_grp_ts.rename(columns={'obs_map':'OBSNME', 'value':'OBSVAL'}, inplace=True)            
            obs_grp_ts = obs_grp_ts[['OBGNME', 'OBSVAL','OBSNME','WEIGHT','model']]
            if index == 0:
                PESTobs = obs_grp_ts
            else:
                PESTobs = PESTobs.append(obs_grp_ts)
            # end if 

#        header = ['OBSNME', 'OBSVAL', 'WEIGHT', 'OBGNME', 'model']
#        num_obs = len(obs.keys())
#        PESTobs = pd.DataFrame(columns=header, index=self.obs.keys())
#        PESTobs['OBSNME'] = obs.keys()
#        PESTobs['OBSVAL'] = obs.values()
#        PESTobs['WEIGHT'] = [1.0] * num_obs
#        PESTobs['OBGNME'] = [key] * num_obs
#        PESTobs['model'] = [self.models_ID[0]] * num_obs

   
        if self.csv_copy:
            self._createCSVcopyFromDataFrame(PESTobs, 'PESTobs')
        #end if
        return PESTobs     
   
    def _createCSVcopyFromDataFrame(self, df, name):   
        fname = self.directory + os.path.sep + name + '.csv'
        if os.path.exists(fname):
            print fname + ' file exists already'
        else:
            df.to_csv(os.path.join(self.directory, name + '.csv'), index=False)        
        #end if
   
    def writePESTdict2csv(self):
        with open(os.path.join(self.directory, self.name + '.csv'), 'w') as f:
            w = csv.DictWriter(f, self.PEST_data.keys())
            w.writeheader()
            w.writerow(self.PEST_data)
#    
#    def readPESTdictFromCsv(self):
#        pass
#
#    def writePESTdict2excel(self):
#        pass
#    
#    def readPESTdictFromExcel(self):
#        pass
#    
    def genParameters(self, method=None):   
        # Generate *name*_parameters.txt containing all models parameters
        #
        # Parameters names and values are taken from 'PARNAME' and 'PARVAL1' in the PESTpar function
        
        # method argument can be either:
        if method:
            if method.lower() in ['csv', 'excel']:
                pass
            else:
                print 'Method not recognized'
                sys.exit(1)
            # end if
        else:
            method = 'dataframe'
        #end if
                
        print 'Generating %s_parameters.txt using %s' %(self.name, method)

        #if self.excel_copy == True:
        #    PESTpar = pd.read_excel('PEST.xlsx', sheetname='PESTpar');
        # end if

        if (method == 'csv'): # & (self.csv_copy == True):
            self.PEST_data['PESTpar'] = pd.read_csv(self.directory + os.path.sep + 'PESTpar.csv')
            self.PEST_data['PESTobs'] = pd.read_csv(self.directory + os.path.sep + 'PESTobs.csv')
        # end if
            
        # Generate *name*_parameters.txt
        self.PEST_data['PESTpar'].to_csv(self.directory + 'parameters.txt', sep='\t', columns=['PARNAME', 'PARVAL1'], index=False)        

    def genPestfiles(self, models_ID=None):     

        # Generate a PEST folder and PEST files for a given chain of models_ID
        # 
        # All necessary PEST inputs are taken from PEST dict containing dicts or dataframes:
        # - control data are read from PESTcon
        # - parameters to calibrate are read from PESTpar, including PARTIED information
        # - parameter groups are read from PESTpgp
        # - observations are read from PESTobs
        # - observation groups are taken as the unique entries of column OBGNME in PESTobs

        if not models_ID:
            models_ID = ['default']
        
        print('# Generating PEST files, %s #\n' %(datetime.datetime.now()))
        
        PEST_name = 'pest'
        PEST_folder_name = 'PEST'

        for model in models_ID:
            PEST_folder_name += '_' + model
        # end for
        
        PEST_folder = self.directory + os.path.sep #+ PEST_folder_name
        if not os.path.isdir(PEST_folder):
            os.mkdir(PEST_folder)
        # end if

        # Retain only parameters corresponding to current models (uses the column 'models' of Dataframe PESTpar) by assigning 'fixed' to the other ones
        for row in self.PEST_data['PESTpar'].iterrows():
            if row[1]['models'] not in models_ID:
                self.PEST_data['PESTpar'].set_value(row[0], 'PARTRANS', 'fixed')
        #end for
        
        # Retain only observations corresponding to current models (uses the column 'model' of Dataframe PESTobs)
        required = []
        for row in self.PEST_data['PESTobs'].iterrows():
            if row[1]['model'] in models_ID:
                required += [row[0]]
        #end for

        self.PEST_data['PESTobs'] = self.PEST_data['PESTobs'].loc[required]
        
        # Special consideration when PEST is run in prediction mode
        if self.PEST_data['PESTcon']['control_data']['PESTMODE'] == 'prediction':
            if self.PEST_data['PESTobs']['OBGNME'] != 'predict':
                prediction_obs_row = {}                
                prediction_obs_row['OBSNME'] = 'prediction'
                prediction_obs_row['OBSVAL'] = 0.0
                prediction_obs_row['WEIGHT'] = 1.0
                prediction_obs_row['OBGNME'] = 'predict'
                prediction_obs_row['model'] = 'any'
                # Add observation to 
                self.PEST_data['PESTobs'].append(prediction_obs_row)
                print('  ** prediction mode was detected\n' +
                    '  --> an observation ''prediction'' was automatically added that belongs to a group called ''predict'', conformally with what PEST requires\n' +
                    '  --> make sure to write the computed prediction at the last line of the observation file which must be writen after your model run\n')
            #end if
        #end if
        
        # Find observation groups from unique observation group names
        PESTobsgp = self.PEST_data['PESTobs']['OBGNME'].unique()

        # Definition of prior information / rules
        
#        """        
#        % PROTOCOL
#        % No more than 300 chars per line, use & + space to continue no next line
#        % Each item separte by at least one space from its neighbors
#        % Start with label <8 chars, case insensitive like all other char
#        % variables.
#        %
#        % To the left of the "=" sign: one or more combinations of
#        % PIFAC * PARNAME|log(PARNAME) + ...
#        % All PARNAME use must be adjustable parameters.
#        % Each parameter can be referenced only once in PRIOR equation.
#        % The parameter factor must be supplied.
#        % log(PARNAME) is necessary if PARNAME is log transformed.
#        % logbase is 10.
#        %
#        % To the right side of the "=" sign are two real variables PIVAL and WEIGHT.
#        % PIVAL is the value that the left side of the the equation aimes to achieve
#        % WEIGHT is the weight of this prior information rule/article.
#        % WEIGHT should preferably be inversly proportional to the standard
#        % deviation of PIVAL, but must always be >=0
#        %
#        % No two prior information articles must say the same thing.
#        %
#        
#        % Adapt to your model
#        %PRIOR={ % PILBL PIFAC * PARNME + PIFAC * log(PARNME) .... = PIVAL WIEGHT
#        %    };
#        """        
        
        PRIOR={}
        
        # Generate PEST control file
        
        # Command line that PEST excutes to run the model
        # This assumes run.bat exists locally and contains a command like 'python run_and_postprocess_model.py'
        # which further assumes that run_and_postprocess_model.py exists and:
        # - loads the model and updates the parameters based on parameters.txt
        # - runs the model
        # - post-processes model results and writes relevant model outputs for PEST to read
     
        if os.name == 'nt':
            PESTCMD = 'run.bat'
        if os.name == 'posix':
            PESTCMD = './run.sh'
            
        # Model input file
        INFLE = 'parameters.txt'
        
        # Corresponding template file for PEST to know how to write it
        TEMPFLE = 'parameters.tpl'

        OUTFLE = {}
        INSFLE = {}        

        os.chdir(self.directory)

        for obs_gp in self.PEST_data['PESTobs']['OBGNME'].unique():
        #for model in models_ID:
            model = models_ID[0]
            # Find model folder
            #model_folder = self.directory + os.path.sep + 'model_' + model
            
            #model_folder = '.' + os.path.sep + 'model_' + model
            #if not os.path.isdir(model_folder):
            #    sys.exit('Model folder not found for model %s' %(model))
            #end if
            
            # Model observation file (must be created by post-processing of model outputs)

            OUTFLE[obs_gp] = '.' + os.path.sep + 'model_' + model + r'\observations_' + obs_gp + '.txt'
            
            # Corresponding instruction file for PEST to know how to read it
            INSFLE[obs_gp] = 'observations_' + obs_gp + '.ins'
        #end for
        
        
        # PEST control file
        PESTFILE = PEST_name + '.pst'
        
        # Counters
        NPAR = self.PEST_data['PESTpar'].count()['PARNAME']
        NOBS = self.PEST_data['PESTobs'].count()['OBSNME']
        NPARGP = self.PEST_data['PESTpgp'].count()['PARGPNME']
        NPRIOR  = len(PRIOR)
        NOBSGP = len(PESTobsgp)
        NTPLFILE = 1
        NINSFLE = len(INSFLE)
        
        # Open file
        with open(os.path.join(PEST_folder, PESTFILE),'w') as f:
            f.write('pcf\n')
            
            # Control data
            control_data = self.PEST_data['PESTcon']['control_data']
            f.write('* control data\n')
            f.write('%s %s\n' %(control_data['RSTFLE'], control_data['PESTMODE']))
            f.write('%d %d %d %d %d\n' % (NPAR, NOBS, NPARGP, NPRIOR, NOBSGP))
            f.write('%d %d %s %s\n' %(NTPLFILE, NINSFLE, control_data['PRECIS'], 
                                      control_data['DPOINT']))
            f.write('%g %g %g %g %d %d %s %s\n' %(
                    control_data['RLAMBDA1'],
                    control_data['RLAMFAC'],
                    control_data['PHIRATSUF'],
                    control_data['PHIREDLAM'],
                    control_data['NUMLAM'],
                    control_data['JACUPDATE'],
                    control_data['LAMFORGIVE'],
                    control_data['DERFORGIVE']))
            f.write('%g %g %g\n' %(
                    control_data['RELPARMAX'],
                    control_data['FACPARMAX'],
                    control_data['FACORIG']))
            f.write('%g %g %s\n' %(
                    control_data['PHIREDSHW'],
                    control_data['NOPTSWITCH'],
                    control_data['BOUNDSCALE']))
            f.write('%d %g %d %d %g %d\n' %(
                    control_data['NOPTMAX'],
                    control_data['PHIREDSTP'],
                    control_data['NPHISTP'],
                    control_data['NPHINORED'],
                    control_data['RELPARSTP'],
                    control_data['NRELPAR']))
            f.write('%d %d %d %d %s %s %s %s %s\n' %(
                    control_data['ICOV'],
                    control_data['ICOR'],
                    control_data['IEIG'],
                    control_data['IRES'],
                    control_data['JCOSAVE'],
                    control_data['VERBOSEREC'],
                    control_data['JCOSAVEITN'],
                    control_data['REISAVEITN'],
                    control_data['PARSAVEITN']))
            
            # SVD
            svd = self.PEST_data['PESTcon']['singular_value_decomposition']
            f.write('* singular value decomposition\n')
            f.write('%d\n' %(svd['SVDMODE']))
            svd['MAXSING'] = len(self.PEST_data['PESTpar'].index)
            f.write('%d %g\n' %(svd['MAXSING'], 
                                svd['EIGTHRESH']))
            f.write('%d\n' %svd['EIGWRITE'])
            # Parameter groups
            f.write('* parameter groups\n')
            for row in self.PEST_data['PESTpgp'].iterrows():
                [f.write(str(x)+'\t') for x in row[1].tolist()]
                f.write('\n')
            #end for

            # Parameters
            f.write('* parameter data\n')
            for row in self.PEST_data['PESTpar'][['PARNAME', 'PARTRANS', 'PARCHGLIM', 'PARVAL1', 'PARLBND', 'PARUBND', 'PARGP', 'SCALE', 'OFFSET']].iterrows():
                [f.write(str(x) + '\t') for x in row[1].tolist()]
                f.write('\n')
            #end for
            for row in self.PEST_data['PESTpar'][['PARNAME', 'PARTRANS', 'PARTIED']].iterrows():
                if row[1]['PARTRANS'] == 'tied':
                    f.write('%s\t%s' %(row[1]['PARNAME'], row[1]['PARTIED']))
                    f.write('\n')
                #end if
            #end for
            
            # Observation groups
            f.write('* observation groups\n');
            [f.write(str(x) + '\n') for x in PESTobsgp]
            #end for
            
            # Observations
            f.write('* observation data\n')
            for row in self.PEST_data['PESTobs'][['OBSNME', 'OBSVAL', 'WEIGHT', 'OBGNME']].iterrows():
               [f.write(str(x) + '\t') for x in row[1]]
               f.write('\n')         
            #end for
            
            # Command line that pest executes
            f.write('* model command line\n')
            f.write('%s\n' %PESTCMD)
            
            # Model input/output
            f.write('* model input/output\n')
            f.write('%s\t%s\n' %(TEMPFLE, INFLE))
            for obs_gp in self.PEST_data['PESTobs']['OBGNME'].unique():
            #for model in models_ID:
               f.write('%s\t%s\n' %(INSFLE[obs_gp], OUTFLE[obs_gp]))
               #f.write('%s\t%s\n' %(INSFLE[model], OUTFLE[model]))
            #end for
            
            # PRIOR rules
            f.write('* prior information\n')
#            for i in range(1:size(PRIOR,1))
#               # PILBL PIFAC * PARNAME + PIFAC * log(PARNAME) + ... = PIVAL WEIGHT'\n');
#               f.write('%s' %(PRIOR[i,1]))  # PRIOR LABEL (no spaces in front
#               for j in range(2:size(PRIOR,2)): # Rest of prior info
#                   if ~isempty(PRIOR(i,j))
#                       if isnum( PRIOR{i,j}), fprintf(fid,' %g',PRIOR{i,j}); end
#                       if ischar(PRIOR{i,j}), fprintf(fid,' %s',PRIOR{i,j}); end
#                   end
#               end
#               f.write('\n')
#            end
            
            # Predictive analysis
            if self.predictive_analysis:
                predictive_analysis = self.PEST_data['PESTcon']['predictive_analysis']
                f.write('* predictive analysis\n')
                f.write('%d\n' %(predictive_analysis['NPREDMAXMIN']))
                f.write('%g %g %g\n' %(predictive_analysis['PD0'], 
                                       predictive_analysis['PD1'], 
                                       predictive_analysis['PD2']))
                f.write('%g %g %g %g %d\n' %(predictive_analysis['ABSPREDLAM'],
                                             predictive_analysis['RELPREDLAM'],
                                             predictive_analysis['INITSCHFAC'],
                                             predictive_analysis['MULSCHFAC'],
                                             predictive_analysis['NSEARCH']))
                f.write('%g %g\n' %(predictive_analysis['ABSPREDSWH'],
                                    predictive_analysis['RELPREDSWH']))
                f.write('%d %g %g %d\n' %(predictive_analysis['NPREDNORED'],
                                          predictive_analysis['ABSPREDSTP'],
                                          predictive_analysis['RELPREDSTP'],
                                          predictive_analysis['NPREDSTP']))
        # end with
        
        
        # Generate initial parameters file
        
#        writetable(table(PESTpar.PARNAME,PESTpar.PARVAL1,'VariableNames',{'PARNAME','PARVAL'}),INFLE{1},'Delimiter','\t');
        
        # Generate PEST template file
        
        with open(os.path.join(PEST_folder, TEMPFLE), 'w') as f:
            f.write('ptf #\n')
            f.write('PARNAME\tPARVAL\n')
            for row in self.PEST_data['PESTpar'][['PARNAME']].iterrows():
                f.write('%s\t#%-15s#\n' %(row[1]['PARNAME'], 
                                              row[1]['PARNAME']))
            #end for
         #end with
        
        # Generate PEST instruction file
        
        for obs_gp in self.PEST_data['PESTobs']['OBGNME'].unique():
        #for model in models_ID: # n in range(1:length(INSFLE)):
            #used_observations = []
            
            #for row in self.PEST_data['PESTobs'].iterrows():
            #    if row[1]['model'] == model:
            #        used_observations += [row[0]]
            #end for
            
            #PESTobs_filtered = self.PEST_data['PESTobs'].loc[used_observations]
            obs_pd = self.PEST_data['PESTobs']
            PESTobs_filtered = obs_pd[obs_pd['OBGNME'] == obs_gp]           

            #with open(PEST_folder + os.path.sep + INSFLE[model],'w') as f:
            with open(os.path.join(PEST_folder, INSFLE[obs_gp]),'w') as f:
                f.write('pif %%\n')
                for row in PESTobs_filtered.iterrows(): #i in range(1:size(PESTobs_n,1))
                    f.write('l1 !%s!\n' %(row[1]['OBSNME']))
                #end for
            #end with
        #end for
        
        # Generate a parameter uncertainty file using bound values
        
        # Calculate parameter standard deviation assuming normal distribution and that lower and upper bounds are 95% intervals
        STD = []        
        for row in self.PEST_data['PESTpar'].iterrows(): #
            # If we assume that the upper and lower bounds specified represent
            # the lower and upper quartiles (i.e. the 25th and 75th percentiles)
            # of the likely parameter distribution,
            # which are + and - 2*sigma, then we can calculate the standard
            # deviation as 1/4 ( upper - lower) = 1/4 ( mean + 2sigma - (mean - 2sigma)) 
            # = 1/4 (4sigma) = sigma            
            if row[1]['PARTRANS'] == 'log':
                #log_trans_params += [row[0]]
                STD += [0.25 * (np.log10(row[1]['PARUBND']) - np.log10(row[1]['PARLBND']))]
            else:
                STD += [0.25 * (row[1]['PARUBND'] - row[1]['PARLBND'])]
            #end if
        #end for                
        self.PEST_data['PESTpar']['STD'] = STD
        
        UNCERTAINTYFILE = PEST_name + '.unc'
        with open(os.path.join(PEST_folder, UNCERTAINTYFILE), 'w') as f:
            f.write('# Parameter uncertainty file\n')
            f.write('# for filling C(k) \n')
            f.write('\n')
#            if pilot_points:
#                for zone in zones:
#                    f.write('START COVARIANCE_MATRIX \n')
#                    f.write('file "points{}.mat" \n'.format(zones))
#                    f.write('variance_multiplier 1.0 \n')
#                    f.write('#first_parameter K1pp1 \n')
#                    f.write('#last_parameter K1pp49 \n')
#                    f.write('END COVARIANCE_MATRIX \n \n')
#                # end for
#            # end if
#            f.write('\n')
            
            # Needs to be for every other parameter
            f.write('START STANDARD_DEVIATION\n')
            for row in self.PEST_data['PESTpar'].iterrows(): #i in range(1:size(PESTpar,1)):
                if row[1]['PARTRANS'] != 'fixed':
                    f.write('%s\t%g\n' %(row[1]['PARNAME'],
                                         row[1]['STD']))
                #end if
            #end for
            f.write('END STANDARD_DEVIATION\n')
        #end with
        
        
#        # Copy necessary files
#        copyfile('run.m',[PEST_folder '\run.m']);
#        copyfile('run.bat',[PEST_folder '\run.bat']);
#        copyfile('runpwtadj1.bat',[PEST_folder '\runpwtadj1.bat']);
#        copyfile('search_replace.py',[PEST_folder '\search_replace.py']);
#        
#
#        # Add regularisation with preferred value
#        system(['"C:\Workspace\bres0010\PEST_Course\USB Stick\pest\addreg1.exe" ' PEST_folder '\' PESTFILE ' ' PEST_folder '\' PEST_name '_reg.pst']);
#        delete([PEST_folder '\' PESTFILE]);
#        movefile([PEST_folder '\' PEST_name '_reg.pst'],[PEST_folder '\' PESTFILE]);
#
#        
#        # Generate pest run files
#        pestbin_folder = 'C:\\Workspace\\bres0010\\PEST_Course\\USB Stick\\pest';
#        beopestbin_folder = 'C:\\Workspace\\bres0010\\PEST_Course\\USB Stick\\beopest';
#        
#        # runpestchek.bat
#        file_text = '@echo off\n' +
#                        '\n' +
#                        ':: running pestchek\n' +
#                        '"' + pestbin_folder + '\\pestchek.exe" pest_pwtadj1.pst\n' +
#                        '\n' +
#                        ':: checking template file\n' +
#                        '"' + pestbin_folder + '\\tempchek.exe" ' + TEMPFLE[1] + '\n' +
#                        '\n' +
#                        ':: checking instruction file\n'
#        for n in range(1:length(INSFLE)):
#           file_text = file_text + '"' + pestbin_folder + '\\inschek.exe" ' + INSFLE[n] '\n'
#        #end for
#        file_name = 'runpestchek.bat'
#        with open(PEST_folder + os.path.sep + file_name,'w') as f:
#            f.write(file_text)
#        #end with
#        
#        # runpest.bat
#        file_text = ['@echo off\n' ...
#                        '\n' ...
#                        '::run pest chek programs\n' ...
#                        'call runpestchek.bat\n' ...
#                        '\n' ...
#                        '::run pest\n' ...
#                        '"' pestbin_folder '\\pest.exe" pest.pst\n'];
#        file_name = 'runpest.bat';
#        fid=fopen([PEST_folder '\' file_name],'wt');
#        fprintf(fid,file_text);
#        fclose(fid);
#        
#        # restartpest.bat
#        file_text = ['@echo off\n' ...
#                        '\n' ...
#                        '::run pest chek programs\n' ...
#                        '::call runpestchek.bat\n' ...
#                        '\n' ...
#                        '::run pest\n' ...
#                        '"' pestbin_folder '\\pest.exe" pest.pst /j\n'];
#        file_name = 'restartpest.bat';
#        fid=fopen([PEST_folder '\' file_name],'wt');
#        fprintf(fid,file_text);
#        fclose(fid);
#        
#        # runbeopest_master.bat
#        file_text = ['@echo off\n' ...
#                        '\n' ...
#                        '::run beopest master\n' ...
#                        '"' beopestbin_folder '\\beopest64.exe" pest_pwtadj1.pst /p1 /H :4004\n'];
#        file_name = 'runbeopest_master.bat';
#        fid=fopen([PEST_folder '\' file_name],'wt');
#        fprintf(fid,file_text);
#        fclose(fid);
#        
#        # runbeopest_slave.bat
#        file_text = ['@echo off\n' ...
#                        '\n' ...
#                        '::run beopest slave\n' ...
#                        '"' beopestbin_folder '\\beopest64.exe" pest_pwtadj1.pst /p1 /H ncgrt70862.isd.ad.flinders.edu.au:4004\n'];
#        file_name = 'runbeopest_slave.bat';
#        fid=fopen([PEST_folder '\' file_name],'wt');
#        fprintf(fid,file_text);
#        fclose(fid);
#        
#        # restartbeopest_master.bat
#        file_text = ['@echo off\n' ...
#                        '\n' ...
#                        '::restart beopest master\n' ...
#                        '"' beopestbin_folder '\\beopest64.exe" pest_pwtadj1.pst /H /s :4004\n'];
#        file_name = 'restartbeopest_master.bat';
#        fid=fopen([PEST_folder '\' file_name],'wt');
#        fprintf(fid,file_text);
#        fclose(fid);
#        
#        # restartbeopest_slave.bat
#        file_text = ['@echo off\n' ...
#                        '\n' ...
#                        '::restart beopest slave\n' ...
#                        '"' beopestbin_folder '\\beopest64.exe" pest_pwtadj1.pst /H /s ncgrt70862.isd.ad.flinders.edu.au:4004\n'];
#        file_name = 'restartbeopest_slave.bat';
#        fid=fopen([PEST_folder '\' file_name],'wt');
#        fprintf(fid,file_text);
#        fclose(fid);
#        
#        
        # Bye-bye message
        
        print('\nPEST files generated\n')
        print(' %s\n' %PESTFILE)
#        # for i=1:length(TEMPFLE), fprintf(' %s\n',TEMPFLE{i}); end;
#        # for i=1:length(INFLE)  , fprintf(' %s\n',INFLE{i});   end;
#        # for i=1:length(INSFLE) , fprintf(' %s\n',INSFLE{i});  end;
#        # fprintf(' %s\n',UNCERTAINTYFILE);
#        # fprintf(' %s\n',[PEST_folder '\run.m']);
#        # fprintf(' %s\n',[PEST_folder '\run.bat']);
#        # fprintf(' %s\n',[PEST_folder '\runpestchek.bat']);
#        # fprintf(' %s\n',[PEST_folder '\runpest.bat']);
#        # fprintf('\n ***** TO RUN PEST --> runpest.bat *****\n');
        print('\nPEST files generation completed!\n')
#        print('\n ******* INSTRUCTIONS *******\n')
#        print('\n 1. runpest.bat --> normally with NOPTMAX=0 to estimate group contributions to objective function\n')
#        print('\n 2. runpwtadj1.bat --> creates pest_pwtadj1.pst with adjusted weights so that all group contributions equal 1\n')

    def updateparameterswithpestbestpar(pestparfile):
        # Generate parameters.txt containing all models parameters
        #
        # Parameters names and values are taken from 'PARNAME' and 'PARVAL1' in PESTpar
        
        print('Generating parameters.txt\n')
        
        # Check existence of the PEST .par file and read it
        
        if not os.path.exists(pestparfile):
            sys.exit("Can't find PEST .par file <<%s>>" %pestparfile)
        #end if

        with open(pestparfile, 'r') as f:
            text = f.readlines()
        #end with
            
#        PESTbestpar = textscan(fileID,'%s %f %f %f\n','HeaderLines',1,'Delimiter',' ','MultipleDelimsAsOne',1);
#        if isnan(PESTbestpar{1,2})
#            fclose(fileID);
#            fileID = fopen(pestparfile);
#            PESTbestpar = textscan(fileID,'%s %f %f %f','HeaderLines',1,'Delimiter','\t','MultipleDelimsAsOne',1);
#        end
        
        
        # Assign PEST best parameter value
        
        self.PEST_data['PESTpar']['PARVAL1'] = PESTbestpar
        
        # Generate *name*_parameters.txt
        
        self.PEST_data['PESTpar'].to_csv(self.directory + 'best_parameters.txt', sep='\t', columns=['PARNAME', 'PARVAL1'], index=False)

#if __name__ == '__main__':
#
#    params = {'p1':12, 'p2':34, 'p3':32}
#    
#    obs = {'ob1':2.0, 'ob2':3.5}
#
#    test = PESTInterface(directory=r"C:\Workspace\part0075\MDB modelling\testbox\PESTtest\\", params=params, obs=obs, csv_copy=True)    
#    
#    test.genParameters(method='csv')
#    
#    
#    test.genPESTpgp()
#    
#    test.genPestfiles(models_ID=['default'])
    