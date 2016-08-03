"""
Interface for PEST software to generate necessary files
Based on code of Etienne Bresciani (applied for mflab in Matlab) translated to 
python and modified as well to work with the GWModelBuilder class.
"""

import os
import sys
import datetime

import pandas as pd

class PESTInterface(object):

    def __init__(self, name=None, directory=None, csv_copy=False, excel_copy=False, param_name_values=None):
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
            
        self.PEST_data['PESTcon'] = self.PESTcon()
        self.PEST_data['PESTpar'] = self.PESTpar()
        self.PEST_data['PESTpgp'] = self.PESTpgp()
        self.PEST_data['PESTobs'] = self.PESTobs()
        
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

    def PESTpar(self):
        header = ['PARNAME', 'PARTRANS', 'PARCHGLIM', 'PARVAL1', 'PARLBND', 'PARUBND', 'PARGP', 'SCALE', 'OFFSET', 'PARTIED', 'models', 'unit', 'comment']
        # For column 'PARNAME' the maximum length is 12 characters
        PARTRANSoptions = ['log', 'fixed', 'tied']
        PARCHGLIMoptions = ['factor', 'relative']
        
        self.PEST_data['PESTpar'] = {}
        self.PEST_data['PESTpar']['header'] = header
        self.PEST_data['PESTpar']['PARTRANSoptions'] = PARTRANSoptions
        self.PEST_data['PESTpar']['PARCHGLIMoptions'] = PARCHGLIMoptions
    
    def PESTpgp(self):
        header = ['PARGPNME', 'INCTYP', 'DERINC', 'DERINCLB', 'FORCEN', 'DERINCMUL', 'DERMTHD']
        INCTYPEdefault = 'relative'
        DERINCdefault = 0.01
        FORCENdefault = 'switch'
        DERMTHDdefault = 'parabolic'        
        
        self.PEST_data['PESTpgp'] = {}
        self.PEST_data['PESTpgp']['header'] = header
        self.PEST_data['PESTpgp']['INTYPEdefault'] = INTYPEdefault
        self.PEST_data['PESTpgp']['DERINCdefault'] =  DERINCdefault
        self.PEST_data['PESTpgp']['FORCENdefault'] =  FORCENdefault
        self.PEST_data['PESTpgp']['DERMTHDdefault'] =  DERMTHDdefault       
        
    def PESTobs(self):
        header = ['OBSNME', 'OBSVAL', 'WEIGHT', 'OBGNME', 'model']
        
        self.PEST_data['PESTobs'] = {}
        self.PEST_data['PESTobs']['header'] = header
    
    def writePESTdict2csv(self):
        pass
    
    def readPESTdictFromCsv(self):
        pass

    def writePESTdict2excel(self):
        pass
    
    def readPESTdictFromExcel(self):
        pass
    
    def genparameters(self, method):   
        # Generate *name*_parameters.txt containing all models parameters
        #
        # Parameters names and values are taken from 'PARNAME' and 'PARVAL1' in the PESTpar function
        
        # method argument can be either:
        if method.lower() in ['class', 'csv', 'excel']:
            pass
        else:
            print 'Method not recognized'
            sys.exit(1)
        # end if
        
        print 'Generating %s_parameters.txt\n using %s' %(self.name, method)
        
        dictfile = directory + self.name + 'PEST.dict'
        if not os.path.isfile(dictfile):
            print "Can't find dict file: %s" %(dictfile)
            sys.exit(1) 
        #end if
    
        if (method == 'class'):
            PESTpar = pd.Dataframe(self.param_name_values)
        # end if


        #if self.excel_copy == True:
        #    PESTpar = pd.read_excel('PEST.xlsx', sheetname='PESTpar');
        # end if

        if (method == 'csv') & (self.csv_copy == True):
            PESTpar = pd.read_csv('PESTpar.csv')
        # end if
            
        # Generate *name*_parameters.txt
        
        PESTvar.to_csv(directory + 'parameters.txt', sep='\t', columns=['PARNAME', 'PARVAL1'])        

    def genpestfiles(self, models_ID):     

        # Generate a PEST folder and PEST files for a given chain of models_ID
        # 
        # All necessary PEST inputs are taken from PEST dict containing dicts or dataframes:
        # - control data are read from PESTcon
        # - parameters to calibrate are read from PESTpar, including PARTIED information
        # - parameter groups are read from PESTpgp
        # - observations are read from PESTobs
        # - observation groups are taken as the unique entries of column OBGNME in PESTobs
        
        print('# Generating PEST files, %s #\n' %(datetime.date()))
        
        
        # Check existence of the dict file
        
        dictfile = directory + self.name + 'PEST.dict'
        if not os.path.isfile(dictfile):
            print "Can't find dict file: %s" %(dictfile)
            sys.exit(1) 
        #end if
        
        
        # Create a folder for this chain of models
        
        current_folder_path = os.getcwd()
        os.chdir('..')
        root_folder = os.getcwd()
        os.chdir(current_folder_path)
        
        PEST_name = 'pest'
        PEST_folder_name = []

        for n in range(1,models_ID+1,2):
            PEST_folder_name = [PEST_folder_name + '_' + str(n)]
        #end for
        
        PEST_folder = root_folder + os.path.sep + PEST_folder_name
        if ~os.path.exist(PEST_folder):
            os.mkdir(PEST_folder)
        #end if
        
        
        # Read relevant worksheets
        
        PESTcon = readtable(excelfile,'filetype','spreadsheet','sheet','PESTcon','readrownames',true)
        PESTpar = readtable(excelfile,'filetype','spreadsheet','sheet','PESTpar')
        PESTpgp = readtable(excelfile,'filetype','spreadsheet','sheet','PESTpgp')
        PESTobs = readtable(excelfile,'filetype','spreadsheet','sheet','PESTobs')
        
        # Retain only parameters corresponding to current models (uses the column 'models' of sheet PESTpar) by assigning 'fixed' to the other ones
        
        used_parameters = false(height(PESTpar),1);
        for i in range(1, len(PESTpar)):
            used_parameters(i) = is_used(PESTpar.models[i],models_ID)
        #end for
            
        PESTpar.PARTRANS(~used_parameters) = 'fixed'
        
        
        # Retain only observations corresponding to current models (uses the column 'model' of sheet PESTobs)
        
        used_observations = false(height(PESTobs),1)
        for i in range(1, len(PESTobs)):
            used_observations(i) = is_used(PESTobs.model[i],models_ID)
        #end for
        PESTobs = PESTobs[used_observations,:]
        
        
        # Special consideration when PEST in run in prediction mode
        
        if PESTcon['control_data']['PESTMODE'] == 'prediction':
            if PESTobs.OBGNME(:) != 'predict':
                prediction_obs_row.OBSNME = 'prediction'
                prediction_obs_row.OBSVAL = 0.0
                prediction_obs_row.WEIGHT = 1.0
                prediction_obs_row.OBGNME = 'predict'
                prediction_obs_row.model = 'any'
                PESTobs = [PESTobs;struct2table(prediction_obs_row)]
                fprintf(['  ** prediction mode was detected\n'...
                    '  --> an observation ''prediction'' was automatically added that belongs to a group called ''predict'', conformally with what PEST requires\n'...
                    '  --> make sure to write the computed prediction at the last line of the observation file which must be writen after your model run\n']);
            #end if
        #end if
        
        
        # Find observation groups from unique observation group names
        
        PESTobsgp = unique(PESTobs.OBGNME)
        
        
        # Definition of prior information / rules
        
        """        
        % PROTOCOL
        % No more than 300 chars per line, use & + space to continue no next line
        % Each item separte by at least one space from its neighbors
        % Start with label <8 chars, case insensitive like all other char
        % variables.
        %
        % To the left of the "=" sign: one or more combinations of
        % PIFAC * PARNAME|log(PARNAME) + ...
        % All PARNAME use must be adjustable parameters.
        % Each parameter can be referenced only once in PRIOR equation.
        % The parameter factor must be supplied.
        % log(PARNAME) is necessary if PARNAME is log transformed.
        % logbase is 10.
        %
        % To the right side of the "=" sign are two real variables PIVAL and WEIGHT.
        % PIVAL is the value that the left side of the the equation aimes to achieve
        % WEIGHT is the weight of this prior information rule/article.
        % WEIGHT should preferably be inversly proportional to the standard
        % deviation of PIVAL, but must always be >=0
        %
        % No two prior information articles must say the same thing.
        %
        
        % Adapt to your model
        %PRIOR={ % PILBL PIFAC * PARNME + PIFAC * log(PARNME) .... = PIVAL WIEGHT
        %    };
        """        
        
        PRIOR={}
        
        
        # Generate PEST control file
        
        # Command line that PEST excutes to run the model
        # This assumes run.bat exists locally and contains a command like 'matlab -wait -nosplash -nodesktop -r run -logfile run.log -minimize'
        # which further assumes that run.m exists and:
        # - set the paths for mfLab
        # - start mf_setup
        # - start a post-processing script that writes model outputs for PEST to read
        PESTCMD = 'run.bat'
        
        # Model input file
        INFLE = '..\common_basis\parameters.txt'
        
        # Corresponding template file for PEST to know how to write it
        TEMPFLE = 'parameters.tpl'
        
        for n in range(1:size(models_ID,2)):
            # Find model folder
            model_folder = ls(root_folder + os.path.sep + models_ID[n] + '_*')
            if ~os.path.exists(model_folder):
                sys.exit('Model folder not found for model %s' %(models_ID[n]))
            #elif size(model_folder,1)>1
            #    error('Several model folders found for model %s', models_ID{n});
            #end
            
            # Model observation file (must be create by post-processing of model outputs)
            OUTFLE[n] = '..' + os.path.sep + model_folder + r'\observations_' + models_ID[n] + '.txt'
            
            # Corresponding instruction file for PEST to know how to read it
            INSFLE[n] = 'observations_' + models_ID[n] + '.ins'
        #end for
        
        
        # PEST control file
        PESTFILE = PEST_name + '.pst'
        
        # Counters
        NPAR    = size(PESTpar,1)
        NOBS    = size(PESTobs,1)
        NPARGP  = size(PESTpgp,1)
        NPRIOR  = size(PRIOR,1)
        NOBSGP  = size(PESTobsgp,1)
        NTPLFILE= len(TEMPFLE)
        NINSFLE = len(INSFLE)
        
        # Open file
        with open(PEST_folder + os.path.sep + PESTFILE,'w') as f:
            f.write('pcf\n')
            
            # Control data
            f.write('* control data\n')
            f.write('%s %s\n' %(PESTcon.VALUE['RSTFLE'], PESTcon.VALUE['PESTMODE']))
            f.write('%d %d %d %d %d\n',NPAR,NOBS,NPARGP,NPRIOR,NOBSGP);
            f.write('%d %d %s %s\n',NTPLFILE,NINSFLE,PESTcon.VALUE{'PRECIS'},PESTcon.VALUE{'DPOINT'});
            f.write('%g %g %g %g %d %d %s %s\n',...
                str2double(PESTcon.VALUE{'RLAMBDA1'}),...
                str2double(PESTcon.VALUE{'RLAMFAC'}),...
                str2double(PESTcon.VALUE{'PHIRATSUF'}),...
                str2double(PESTcon.VALUE{'PHIREDLAM'}),...
                str2double(PESTcon.VALUE{'NUMLAM'}),...
                str2double(PESTcon.VALUE{'JACUPDATE'}),...
                PESTcon.VALUE{'LAMFORGIVE'},...
                PESTcon.VALUE{'DERFORGIVE'});
            f.write('%g %g %g\n',...
                str2double(PESTcon.VALUE{'RELPARMAX'}),...
                str2double(PESTcon.VALUE{'FACPARMAX'}),...
                str2double(PESTcon.VALUE{'FACORIG'}));
            f.write('%g %g %s\n',...
                str2double(PESTcon.VALUE{'PHIREDSHW'}),...
                str2double(PESTcon.VALUE{'NOPTSWITCH'}),...
                PESTcon.VALUE{'BOUNDSCALE'});
            f.write('%d %g %d %d %g %d\n',...
                str2double(PESTcon.VALUE{'NOPTMAX'}),...
                str2double(PESTcon.VALUE{'PHIREDSTP'}),...
                str2double(PESTcon.VALUE{'NPHISTP'}),...
                str2double(PESTcon.VALUE{'NPHINORED'}),...
                str2double(PESTcon.VALUE{'RELPARSTP'}),...
                str2double(PESTcon.VALUE{'NRELPAR'}));
            f.write('%d %d %d\n',...
                str2double(PESTcon.VALUE{'ICOV'}),...
                str2double(PESTcon.VALUE{'ICOR'}),...
                str2double(PESTcon.VALUE{'IEIG'}));
            
            # SVD
            f.write('* singular value decomposition\n');
            f.write('%d\n',str2double(PESTcon.VALUE{'SVDMODE'}));
            f.write('%d %g\n',str2double(PESTcon.VALUE{'MAXSING'}),str2double(PESTcon.VALUE{'EIGTHRESH'}));
            f.write('%d\n',str2double(PESTcon.VALUE{'EIGWRITE'}));
            
            # Parameter groups
            f.write('* parameter groups\n');
            for i in range(1:size(PESTpgp,1)):
               f.write('%s\t%s\t%g\t%g\t%s\t%g\t%s\n',...
                   PESTpgp.PARGPNME{i},PESTpgp.INCTYP{i},PESTpgp.DERINC(i),PESTpgp.DERINCLB(i),PESTpgp.FORCEN{i},...
                   PESTpgp.DERINCMUL(i),PESTpgp.DERMTHD{i});
            #end for
            
            # Parameters
            f.write('* parameter data\n');
            for i in range(1:size(PESTpar,1)):
                f.write('%s\t%s\t%s\t%g\t%g\t%g\t%s\t%g\t%g\n',...
                   PESTpar.PARNAME{i},PESTpar.PARTRANS{i},PESTpar.PARCHGLIM{i},PESTpar.PARVAL1(i),PESTpar.PARLBND(i),...
                   PESTpar.PARUBND(i),PESTpar.PARGP{i},PESTpar.SCALE(i),PESTpar.OFFSET(i));
            #end for
            for i in range(1:size(PESTpar,1)):
                if strcmpi(PESTpar{i,2},'TIED')
                    f.write('%s\t%s\n',PESTpar.PARNAME{i},PESTpar.PARTIED{i});
                #end if
            #end for
            
            # Observation groups
            f.write('* observation groups\n');
            for i in range(1:size(PESTobsgp,1))
               f.write('%s\n' %(PESTobsgp[i]))
            #end for
            
            # Observations
            f.write('* observation data\n')
            for i in range(1:size(PESTobs,1)):
               f.write('%s\t%g\t%g\t%s\n',...
                   PESTobs.OBSNME{i},PESTobs.OBSVAL(i),PESTobs.WEIGHT(i),PESTobs.OBGNME{i});
            #end for
            
            # Command line that pest executes
            f.write('* model command line\n');
            f.write('%s\n',PESTCMD);
            
            # Model input/output
            f.write('* model input/output\n');
            for i in range(1:length(TEMPFLE)):
               f.write('%s\t%s\n',TEMPFLE{i},INFLE{i})
            #end for
            for i in range(1:length(INSFLE)):
               f.write('%s\t%s\n',INSFLE{i},OUTFLE{i});
            #end for
            
            # PRIOR rules
            f.write('* prior information\n')
            for i in range(1:size(PRIOR,1))
               # PILBL PIFAC * PARNAME + PIFAC * log(PARNAME) + ... = PIVAL WEIGHT'\n');
               f.write('%s' %(PRIOR[i,1]))  # PRIOR LABEL (no spaces in front
               for j in range(2:size(PRIOR,2)): # Rest of prior info
                   if ~isempty(PRIOR(i,j))
                       if isnum( PRIOR{i,j}), fprintf(fid,' %g',PRIOR{i,j}); end
                       if ischar(PRIOR{i,j}), fprintf(fid,' %s',PRIOR{i,j}); end
                   end
               end
               f.write('\n')
            end
            
            # Predictive analysis
            f.write('* predictive analysis\n')
            f.write('%d\n' %(PESTcon.VALUE['NPREDMAXMIN']))
            f.write('%g %g %g\n' %(PESTcon.VALUE['PD0'], 
                                   PESTcon.VALUE['PD1'], 
                                   PESTcon.VALUE['PD2']))
            f.write('%g %g %g %g %d\n' %(PESTcon.VALUE['ABSPREDLAM'],
                                         PESTcon.VALUE['RELPREDLAM'],
                                         PESTcon.VALUE['INITSCHFAC'],
                                         PESTcon.VALUE['MULSCHFAC'],
                                         PESTcon.VALUE['NSEARCH']))
            f.write('%g %g\n' %(PESTcon.VALUE['ABSPREDSWH'],
                                PESTcon.VALUE['RELPREDSWH']))
            f.write('%d %g %g %d\n' %(PESTcon.VALUE['NPREDNORED'],
                                      PESTcon.VALUE['ABSPREDSTP'],
                                      PESTcon.VALUE['RELPREDSTP'],
                                      PESTcon.VALUE['NPREDSTP']))
        # end with
        
        
        # Generate initial parameters file
        
        writetable(table(PESTpar.PARNAME,PESTpar.PARVAL1,'VariableNames',{'PARNAME','PARVAL'}),INFLE{1},'Delimiter','\t');
        
        
        # Generate PEST template file
        
        with open(PEST_folder + os.path.sep + TEMPFLE[1],'w') as f:
            f.write('ptf #\n')
            f.write('PARNAME\tPARVAL\n')
            for i in range(1:size(PESTpar,1)):
                fprintf(fid,'%s\t#%-15s#\n' %(PESTpar.PARNAME[i], 
                                              PESTpar.PARNAME[i])
            #end for
         #end with
        
        # Generate PEST instruction file
        
        for n in range(1:length(INSFLE)):
            used_observations_n = false(height(PESTobs),1)
            for i in range(1:height(PESTobs)):
                used_observations_n(i) = is_used(PESTobs.model{i},models_ID(n))
            #end for
            PESTobs_n = PESTobs(used_observations_n,:);
            with open(PEST_folder + os.path.sep + INSFLE[n],'w') as f:
                f.write('pif %%\n')
                for i in range(1:size(PESTobs_n,1))
                    f.write('l1 !%s!\n' %(PESTobs_n.OBSNME[i]))
                #end for
            #end with
        #end for
        
        # Generate a parameter uncertainty file using bound values
        
        ## Calculate parameter standard deviation assuming normal distribution and that lower and upper bounds are 95% intervals
        PESTpar.STD = 0.25 * (PESTpar.PARUBND - PESTpar.PARLBND);
        log_trans_params = strcmp(PESTpar.PARTRANS,'log');
        PESTpar.STD(log_trans_params) = 0.25 * (log10(PESTpar.PARUBND(log_trans_params)) - log10(PESTpar.PARLBND(log_trans_params)));
        
        UNCERTAINTYFILE = PEST_name + '.unc'
        with open(PEST_folder + os.path.sep + UNCERTAINTYFILE],'w') as f:
            f.write('# Uncertainty file\n');
            f.write('\n');
            f.write('START STANDARD_DEVIATION\n');
            for i in range(1:size(PESTpar,1)):
                if PESTpar.PARTRANS[i] != 'fixed':
                    f.write('%s\t%g\n' %(PESTpar.PARNAME[i]
                                         PESTpar.STD[i]))
                #end if
            #end for
            f.write('END STANDARD_DEVIATION\n')
        #end with
        
        
        # Copy necessary files
        copyfile('run.m',[PEST_folder '\run.m']);
        copyfile('run.bat',[PEST_folder '\run.bat']);
        copyfile('runpwtadj1.bat',[PEST_folder '\runpwtadj1.bat']);
        copyfile('search_replace.py',[PEST_folder '\search_replace.py']);
        

        # Add regularisation with preferred value
        system(['"C:\Workspace\bres0010\PEST_Course\USB Stick\pest\addreg1.exe" ' PEST_folder '\' PESTFILE ' ' PEST_folder '\' PEST_name '_reg.pst']);
        delete([PEST_folder '\' PESTFILE]);
        movefile([PEST_folder '\' PEST_name '_reg.pst'],[PEST_folder '\' PESTFILE]);

        
        # Generate pest run files
        pestbin_folder = 'C:\\Workspace\\bres0010\\PEST_Course\\USB Stick\\pest';
        beopestbin_folder = 'C:\\Workspace\\bres0010\\PEST_Course\\USB Stick\\beopest';
        
        # runpestchek.bat
        file_text = '@echo off\n' +
                        '\n' +
                        ':: running pestchek\n' +
                        '"' + pestbin_folder + '\\pestchek.exe" pest_pwtadj1.pst\n' +
                        '\n' +
                        ':: checking template file\n' +
                        '"' + pestbin_folder + '\\tempchek.exe" ' + TEMPFLE[1] + '\n' +
                        '\n' +
                        ':: checking instruction file\n'
        for n in range(1:length(INSFLE)):
           file_text = file_text + '"' + pestbin_folder + '\\inschek.exe" ' + INSFLE[n] '\n'
        #end for
        file_name = 'runpestchek.bat'
        with open(PEST_folder + os.path.sep + file_name,'w') as f:
            f.write(file_text)
        #end with
        
        # runpest.bat
        file_text = ['@echo off\n' ...
                        '\n' ...
                        '::run pest chek programs\n' ...
                        'call runpestchek.bat\n' ...
                        '\n' ...
                        '::run pest\n' ...
                        '"' pestbin_folder '\\pest.exe" pest.pst\n'];
        file_name = 'runpest.bat';
        fid=fopen([PEST_folder '\' file_name],'wt');
        fprintf(fid,file_text);
        fclose(fid);
        
        # restartpest.bat
        file_text = ['@echo off\n' ...
                        '\n' ...
                        '::run pest chek programs\n' ...
                        '::call runpestchek.bat\n' ...
                        '\n' ...
                        '::run pest\n' ...
                        '"' pestbin_folder '\\pest.exe" pest.pst /j\n'];
        file_name = 'restartpest.bat';
        fid=fopen([PEST_folder '\' file_name],'wt');
        fprintf(fid,file_text);
        fclose(fid);
        
        # runbeopest_master.bat
        file_text = ['@echo off\n' ...
                        '\n' ...
                        '::run beopest master\n' ...
                        '"' beopestbin_folder '\\beopest64.exe" pest_pwtadj1.pst /p1 /H :4004\n'];
        file_name = 'runbeopest_master.bat';
        fid=fopen([PEST_folder '\' file_name],'wt');
        fprintf(fid,file_text);
        fclose(fid);
        
        # runbeopest_slave.bat
        file_text = ['@echo off\n' ...
                        '\n' ...
                        '::run beopest slave\n' ...
                        '"' beopestbin_folder '\\beopest64.exe" pest_pwtadj1.pst /p1 /H ncgrt70862.isd.ad.flinders.edu.au:4004\n'];
        file_name = 'runbeopest_slave.bat';
        fid=fopen([PEST_folder '\' file_name],'wt');
        fprintf(fid,file_text);
        fclose(fid);
        
        # restartbeopest_master.bat
        file_text = ['@echo off\n' ...
                        '\n' ...
                        '::restart beopest master\n' ...
                        '"' beopestbin_folder '\\beopest64.exe" pest_pwtadj1.pst /H /s :4004\n'];
        file_name = 'restartbeopest_master.bat';
        fid=fopen([PEST_folder '\' file_name],'wt');
        fprintf(fid,file_text);
        fclose(fid);
        
        # restartbeopest_slave.bat
        file_text = ['@echo off\n' ...
                        '\n' ...
                        '::restart beopest slave\n' ...
                        '"' beopestbin_folder '\\beopest64.exe" pest_pwtadj1.pst /H /s ncgrt70862.isd.ad.flinders.edu.au:4004\n'];
        file_name = 'restartbeopest_slave.bat';
        fid=fopen([PEST_folder '\' file_name],'wt');
        fprintf(fid,file_text);
        fclose(fid);
        
        
        # Bye-bye message
        
        # fprintf('\nPEST files generated\n');
        # fprintf(' %s\n',PESTFILE);
        # for i=1:length(TEMPFLE), fprintf(' %s\n',TEMPFLE{i}); end;
        # for i=1:length(INFLE)  , fprintf(' %s\n',INFLE{i});   end;
        # for i=1:length(INSFLE) , fprintf(' %s\n',INSFLE{i});  end;
        # fprintf(' %s\n',UNCERTAINTYFILE);
        # fprintf(' %s\n',[PEST_folder '\run.m']);
        # fprintf(' %s\n',[PEST_folder '\run.bat']);
        # fprintf(' %s\n',[PEST_folder '\runpestchek.bat']);
        # fprintf(' %s\n',[PEST_folder '\runpest.bat']);
        # fprintf('\n ***** TO RUN PEST --> runpest.bat *****\n');
        print('\nPEST files generation completed!\n')
        print('\n ******* INSTRUCTIONS *******\n')
        print('\n 1. runpest.bat --> normally with NOPTMAX=0 to estimate group contributions to objective function\n')
        print('\n 2. runpwtadj1.bat --> creates pest_pwtadj1.pst with adjusted weights so that all group contributions equal 1\n')
        
        
        function is_used = is_used(param_models,current_models)
            for k = 1:size(current_models,2)
                if ~isempty(strfind(param_models,current_models{k})) || strcmp(param_models,'any')
                    is_used = true;
                    return;
                end
            end
            is_used = false;

    def updateparameterswithpestbestpar(pestparfile):
        # Generate parameters.txt containing all models parameters
        #
        # Parameters names and values are taken from 'PARNAME' and 'PARVAL1' in PESTpar
        
        print('Generating parameters.txt\n');
        
        # Check existence of the Excel file and read relevant worksheet
        
        excelfile = 'PEST.xlsm';
        if ~exist(excelfile,'file')
            error('Can''t find Excel file <<%s>>',excelfile);
        end
        PESTpar = readtable(excelfile,'filetype','spreadsheet','sheet','PESTpar');
        
        
        # Check existence of the PEST .par file and read it
        
        if ~exist(pestparfile,'file')
            error('Can''t find PEST .par file <<%s>>',pestparfile);
        end
        fileID = fopen(pestparfile);
        PESTbestpar = textscan(fileID,'%s %f %f %f\n','HeaderLines',1,'Delimiter',' ','MultipleDelimsAsOne',1);
        if isnan(PESTbestpar{1,2})
            fclose(fileID);
            fileID = fopen(pestparfile);
            PESTbestpar = textscan(fileID,'%s %f %f %f','HeaderLines',1,'Delimiter','\t','MultipleDelimsAsOne',1);
        end
        fclose(fileID);
        
        
        # Assign PEST best parameter value
        
        PESTpar.PARVAL1 = PESTbestpar{:,2};
        
        
        # Generate *name*_parameters.txt
        
        writetable(table(PESTpar.PARNAME,PESTpar.PARVAL1,'VariableNames',{'PARNAME','PARVAL'}),'..\common_basis\parameters.txt','Delimiter','\t');

    