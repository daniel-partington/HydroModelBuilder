"""
Interface for PEST software to generate necessary files
Based on code of Etienne Bresciani (applied for mflab in Matlab) translated to 
python and modified as well to work with the GWModelBuilder class.
"""

import os
import sys
import datetime
import csv
import pandas as pd

class PESTInterface(object):

    def __init__(self, name=None, directory=None, csv_copy=False, excel_copy=False, params=None, obs=None):
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
            
        self.obs = obs
        
        self.PEST_data['PESTcon'] = {}
        self.PESTcon()
        self.PEST_data['PESTpar'] = self.PESTpar(params=self.params)
        
        #self.PEST_data['PESTpgp'] = {}
        #self.PESTpgp()
        self.PEST_data['PESTobs'] = self.PESTobs(obs=self.obs)
        
    def PESTcon(self):
        control_data = {'RSTFLE': 'restart',
                        'PESTMODE': 'estimation',
                        'PRECIS': 'single',
                        'DPOINT': 'point',
                        'RLAMBDA1': 10,
                        'RLAMFAC': -3,
                        'PHIRATSUF': 0.3,
                        'PHIREDLAM': 1.00E-02,
                        'NUMLAM': -20,
                        'JACUPDATE': 999,
                        'LAMFORGIVE': 'lamforgive',
                        'DERFORGIVE': 'derforgive',
                        'RELPARMAX': 0.5,
                        'FACPARMAX': 5,
                        'FACORIG': 1.00E-04,
                        'PHIREDSHW': 0.02,
                        'NOPTSWITCH': 6,
                        'BOUNDSCALE': 'boundscale',
                        'NOPTMAX':	0,
                        'PHIREDSTP':	0.005,
                        'NPHISTP':	4,
                        'NPHINORED': 4,
                        'RELPARSTP': 0.005,
                        'NRELPAR': 4,
                        'ICOV': 1,
                        'ICOR': 1,
                        'IEIG': 1
                        }        
                        
        singular_value_decomposition = {'SVDMODE': 1,
                                        'MAXSING': 68, # Number of parameters
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
        PESTpar['PARTRANS'] = ['log'] * num_param 
        PESTpar['PARCHGLIM'] = ['factor'] * num_param         
        PESTpar['PARVAL1'] = params.values()
        PESTpar['PARLBND'] = [0] * num_param        
        PESTpar['PARUBND'] = [0] * num_param
        PESTpar['PARGP'] = ['default'] * num_param
        PESTpar['SCALE'] = [1.0] * num_param
        PESTpar['OFFSET'] = [0.0] * num_param
        PESTpar['PARTIED'] = [''] * num_param
        PESTpar['models'] = ['default'] * num_param
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
        
    def PESTobs(self, obs=None):
        header = ['OBSNME', 'OBSVAL', 'WEIGHT', 'OBGNME', 'model']

        num_obs = len(obs.keys())
        PESTobs = pd.DataFrame(columns=header, index=self.obs.keys())
        PESTobs['OBSNME'] = obs.keys()
        PESTobs['OBSVAL'] = obs.values()
        PESTobs['WEIGHT'] = [1.0] * num_obs
        PESTobs['OBGNME'] = ['default'] * num_obs
        PESTobs['model'] = ['default'] * num_obs
   
        if self.csv_copy:
            self._createCSVcopyFromDataFrame(PESTobs, 'PESTobs')
        #end if
        return PESTobs     
   
    def _createCSVcopyFromDataFrame(self, df, name):   
        fname = self.directory + os.path.sep + name + '.csv'
        if os.path.exists(fname):
            print fname + ' file exists already'
        else:
            df.to_csv(self.directory + os.path.sep + name + '.csv', index=False)        
        #end if
   
    def writePESTdict2csv(self):
        with open(self.directory + os.path.sep + self.name +'.csv', 'w') as f:
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
            print row[1]['model']
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
        PESTCMD = 'run.bat'
        
        # Model input file
        INFLE = 'parameters.txt'
        
        # Corresponding template file for PEST to know how to write it
        TEMPFLE = 'parameters.tpl'

        OUTFLE = {}
        INSFLE = {}        
        for model in models_ID:
            # Find model folder
            model_folder = self.directory + os.path.sep + model
            if not os.path.isdir(model_folder):
                sys.exit('Model folder not found for model %s' %(model))
            #end if
            
            # Model observation file (must be create by post-processing of model outputs)
            OUTFLE[model] = '.' + os.path.sep + model + r'\observations_' + model + '.txt'
            
            # Corresponding instruction file for PEST to know how to read it
            INSFLE[model] = 'observations_' + model + '.ins'
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
        with open(PEST_folder + os.path.sep + PESTFILE,'w') as f:
            f.write('pcf\n')
            
            # Control data
            f.write('* control data\n')
            f.write('%s %s\n' %(self.PEST_data['PESTcon']['control_data']['RSTFLE'], self.PEST_data['PESTcon']['control_data']['PESTMODE']))
            f.write('%d %d %d %d %d\n' % (NPAR, NOBS, NPARGP, NPRIOR, NOBSGP))
            f.write('%d %d %s %s\n' %(NTPLFILE, NINSFLE, self.PEST_data['PESTcon']['control_data']['PRECIS'], 
                                      self.PEST_data['PESTcon']['control_data']['DPOINT']))
            f.write('%g %g %g %g %d %d %s %s\n' %(
                    self.PEST_data['PESTcon']['control_data']['RLAMBDA1'],
                    self.PEST_data['PESTcon']['control_data']['RLAMFAC'],
                    self.PEST_data['PESTcon']['control_data']['PHIRATSUF'],
                    self.PEST_data['PESTcon']['control_data']['PHIREDLAM'],
                    self.PEST_data['PESTcon']['control_data']['NUMLAM'],
                    self.PEST_data['PESTcon']['control_data']['JACUPDATE'],
                    self.PEST_data['PESTcon']['control_data']['LAMFORGIVE'],
                    self.PEST_data['PESTcon']['control_data']['DERFORGIVE']))
            f.write('%g %g %g\n' %(
                    self.PEST_data['PESTcon']['control_data']['RELPARMAX'],
                    self.PEST_data['PESTcon']['control_data']['FACPARMAX'],
                    self.PEST_data['PESTcon']['control_data']['FACORIG']))
            f.write('%g %g %s\n' %(
                    self.PEST_data['PESTcon']['control_data']['PHIREDSHW'],
                    self.PEST_data['PESTcon']['control_data']['NOPTSWITCH'],
                    self.PEST_data['PESTcon']['control_data']['BOUNDSCALE']))
            f.write('%d %g %d %d %g %d\n' %(
                    self.PEST_data['PESTcon']['control_data']['NOPTMAX'],
                    self.PEST_data['PESTcon']['control_data']['PHIREDSTP'],
                    self.PEST_data['PESTcon']['control_data']['NPHISTP'],
                    self.PEST_data['PESTcon']['control_data']['NPHINORED'],
                    self.PEST_data['PESTcon']['control_data']['RELPARSTP'],
                    self.PEST_data['PESTcon']['control_data']['NRELPAR']))
            f.write('%d %d %d\n' %(
                    self.PEST_data['PESTcon']['control_data']['ICOV'],
                    self.PEST_data['PESTcon']['control_data']['ICOR'],
                    self.PEST_data['PESTcon']['control_data']['IEIG']))
            
            # SVD
            f.write('* singular value decomposition\n')
            f.write('%d\n' %(self.PEST_data['PESTcon']['singular_value_decomposition']['SVDMODE']))
            f.write('%d %g\n' %(self.PEST_data['PESTcon']['singular_value_decomposition']['MAXSING'], 
                                self.PEST_data['PESTcon']['singular_value_decomposition']['EIGTHRESH']))
            f.write('%d\n' %self.PEST_data['PESTcon']['singular_value_decomposition']['EIGWRITE'])
            
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
            for model in models_ID:
               f.write('%s\t%s\n' %(INSFLE[model], OUTFLE[model]))
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
            f.write('* predictive analysis\n')
            f.write('%d\n' %(self.PEST_data['PESTcon']['predictive_analysis']['NPREDMAXMIN']))
            f.write('%g %g %g\n' %(self.PEST_data['PESTcon']['predictive_analysis']['PD0'], 
                                   self.PEST_data['PESTcon']['predictive_analysis']['PD1'], 
                                   self.PEST_data['PESTcon']['predictive_analysis']['PD2']))
            f.write('%g %g %g %g %d\n' %(self.PEST_data['PESTcon']['predictive_analysis']['ABSPREDLAM'],
                                         self.PEST_data['PESTcon']['predictive_analysis']['RELPREDLAM'],
                                         self.PEST_data['PESTcon']['predictive_analysis']['INITSCHFAC'],
                                         self.PEST_data['PESTcon']['predictive_analysis']['MULSCHFAC'],
                                         self.PEST_data['PESTcon']['predictive_analysis']['NSEARCH']))
            f.write('%g %g\n' %(self.PEST_data['PESTcon']['predictive_analysis']['ABSPREDSWH'],
                                self.PEST_data['PESTcon']['predictive_analysis']['RELPREDSWH']))
            f.write('%d %g %g %d\n' %(self.PEST_data['PESTcon']['predictive_analysis']['NPREDNORED'],
                                      self.PEST_data['PESTcon']['predictive_analysis']['ABSPREDSTP'],
                                      self.PEST_data['PESTcon']['predictive_analysis']['RELPREDSTP'],
                                      self.PEST_data['PESTcon']['predictive_analysis']['NPREDSTP']))
        # end with
        
        
        # Generate initial parameters file
        
#        writetable(table(PESTpar.PARNAME,PESTpar.PARVAL1,'VariableNames',{'PARNAME','PARVAL'}),INFLE{1},'Delimiter','\t');
        
        
        # Generate PEST template file
        
        with open(PEST_folder + os.path.sep + TEMPFLE, 'w') as f:
            f.write('ptf #\n')
            f.write('PARNAME\tPARVAL\n')
            for row in self.PEST_data['PESTpar'][['PARNAME']].iterrows():
                f.write('%s\t#%-15s#\n' %(row[1]['PARNAME'], 
                                              row[1]['PARNAME']))
            #end for
         #end with
        
        # Generate PEST instruction file
        
        for model in models_ID: # n in range(1:length(INSFLE)):
            used_observations = [False] * NOBS
            for i in range(NOBS):
                used_observations(i) = is_used(PESTobs.model{i},models_ID(n))
            #end for
            PESTobs_filtered = self.PEST_data['PESTobs'].loc[used_observations]
            
            with open(PEST_folder + os.path.sep + INSFLE[model],'w') as f:
                f.write('pif %%\n')
                for row in PEST_obs_filtered.iterrows(): #i in range(1:size(PESTobs_n,1))
                    f.write('l1 !%s!\n' %(row[1]['OBSNME']))
                #end for
            #end with
        #end for
        
        # Generate a parameter uncertainty file using bound values
        
        ## Calculate parameter standard deviation assuming normal distribution and that lower and upper bounds are 95% intervals
#        PESTpar.STD = 0.25 * (PESTpar.PARUBND - PESTpar.PARLBND);
#        log_trans_params = strcmp(PESTpar.PARTRANS,'log');
#        PESTpar.STD(log_trans_params) = 0.25 * (log10(PESTpar.PARUBND(log_trans_params)) - log10(PESTpar.PARLBND(log_trans_params)));
#        
#        UNCERTAINTYFILE = PEST_name + '.unc'
#        with open(PEST_folder + os.path.sep + UNCERTAINTYFILE],'w') as f:
#            f.write('# Uncertainty file\n');
#            f.write('\n');
#            f.write('START STANDARD_DEVIATION\n');
#            for i in range(1:size(PESTpar,1)):
#                if PESTpar.PARTRANS[i] != 'fixed':
#                    f.write('%s\t%g\n' %(PESTpar.PARNAME[i]
#                                         PESTpar.STD[i]))
#                #end if
#            #end for
#            f.write('END STANDARD_DEVIATION\n')
#        #end with
#        
#        
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
#        print('\nPEST files generation completed!\n')
#        print('\n ******* INSTRUCTIONS *******\n')
#        print('\n 1. runpest.bat --> normally with NOPTMAX=0 to estimate group contributions to objective function\n')
#        print('\n 2. runpwtadj1.bat --> creates pest_pwtadj1.pst with adjusted weights so that all group contributions equal 1\n')
#        
#
#    def updateparameterswithpestbestpar(pestparfile):
#        # Generate parameters.txt containing all models parameters
#        #
#        # Parameters names and values are taken from 'PARNAME' and 'PARVAL1' in PESTpar
#        
#        print('Generating parameters.txt\n');
#        
#        # Check existence of the Excel file and read relevant worksheet
#        
#        excelfile = 'PEST.xlsm';
#        if ~exist(excelfile,'file')
#            error('Can''t find Excel file <<%s>>',excelfile);
#        end
#        PESTpar = readtable(excelfile,'filetype','spreadsheet','sheet','PESTpar');
#        
#        
#        # Check existence of the PEST .par file and read it
#        
#        if ~exist(pestparfile,'file')
#            error('Can''t find PEST .par file <<%s>>',pestparfile);
#        end
#        fileID = fopen(pestparfile);
#        PESTbestpar = textscan(fileID,'%s %f %f %f\n','HeaderLines',1,'Delimiter',' ','MultipleDelimsAsOne',1);
#        if isnan(PESTbestpar{1,2})
#            fclose(fileID);
#            fileID = fopen(pestparfile);
#            PESTbestpar = textscan(fileID,'%s %f %f %f','HeaderLines',1,'Delimiter','\t','MultipleDelimsAsOne',1);
#        end
#        fclose(fileID);
#        
#        
#        # Assign PEST best parameter value
#        
#        PESTpar.PARVAL1 = PESTbestpar{:,2};
#        
#        
#        # Generate *name*_parameters.txt
#        
#        writetable(table(PESTpar.PARNAME,PESTpar.PARVAL1,'VariableNames',{'PARNAME','PARVAL'}),'..\common_basis\parameters.txt','Delimiter','\t');

if __name__ == '__main__':

    params = {'p1':12, 'p2':34, 'p3':32}
    
    obs = {'ob1':2.0, 'ob2':3.5}

    test = PESTInterface(directory=r"C:\Workspace\part0075\MDB modelling\testbox\PESTtest\\", params=params, obs=obs, csv_copy=True)    
    
    test.genParameters(method='csv')
    
    
    test.genPESTpgp()
    
    test.genPestfiles(models_ID=['default'])
    