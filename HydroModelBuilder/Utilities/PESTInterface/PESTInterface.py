"""
Interface for PEST software to generate necessary files
Based on code of Etienne Bresciani (applied for mflab in Matlab) translated to
python and modified as well to work with the GWModelBuilder class.
"""
import csv
import datetime
import os

import numpy as np
import pandas as pd

from HydroModelBuilder.Utilities.text_writer import write_line, write_multiline


class PESTInterface(object):
    """TODO: Docs"""

    def __init__(self, name=None, directory=None, csv_copy=False,
                 excel_copy=False, params=None, obs=None, obs_grp=None,
                 models_ID=None, predictive_analysis=False):
        self.PEST_data = {}

        self.name = 'default' if name is None else name
        self.directory = '' if directory is None else directory

        self.csv_copy = csv_copy
        self.excel_copy = excel_copy

        self.params = params if params else {}
        self.obs = obs if obs else {}
        self.obs_grp = obs_grp if obs_grp else {}
        self.predictive_analysis = predictive_analysis
        self.models_ID = models_ID

        self.PESTcon()
        self.PEST_data['PESTpar'] = self.PESTpar(params=self.params)
        self.PEST_data['PESTobs'] = self.PESTobs(obs_grp=self.obs_grp)

    def PESTcon(self):
        """Control data for the *.pst file of pest, which takes the form:

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

        PHIREDLAM = ... A suitable value for PHIREDLAM is between 0.01 and 0.05 [float]

        NUMLAM = This integer variable places an upper limit on the number of
                 lambdas that PEST will test during any one iteration. It
                 should normally be set between 5 and 10 (normally closer to 10);
                 however if RLAMBDA1 is set to zero (which is not recommended)
                 it must be set to 1. [integer]

        JACUPDATE = The Broyden Jacobian update procedure is described in
                    section 5.4.2 of Doherty (2015). It provides a mechanism
                    for improving the Jacobian matrix based on model outputs
                    calculated during model runs undertaken for the purpose of
                    testing parameter upgrades calculated using different
                    values of the Marquardt lambda. [integer, optional]

        LAMFORGIVE =

        [DERFORGIVE]
        RELPARMAX FACPARMAX FACORIG [IBOUNDSTICK UPVECBEND] [ABSPARMAX]
        PHIREDSWH [NOPTSWITCH] [SPLITSWH] [DOAUI] [DOSENREUSE] [BOUNDSCALE]
        NOPTMAX PHIREDSTP NPHISTP NPHINORED RELPARSTP NRELPAR [PHISTOPTHRESH] [LASTRUN] [PHIABANDON]
        ICOV ICOR IEIG [IRES] [JCOSAVE] [VERBOSEREC] [JCOSAVEITN] [REISAVEITN] [PARSAVEITN] [PARSAVERUN]
        ...

        *** end file ***

        NOTE: Run pestchek to insure that all control data has been entered appropriately"""

        control_data = {'RSTFLE': 'restart',
                        'PESTMODE': 'estimation',
                        'PRECIS': 'single',
                        'DPOINT': 'point',
                        'RLAMBDA1': 10.,
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
                        'BOUNDSCALE': 'noboundscale',  # 'boundscale',
                        'NOPTMAX': 0,
                        'PHIREDSTP': 0.01,
                        'NPHISTP': 5,
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
                                        'MAXSING': 0,  # Number of parameters
                                        'EIGTHRESH': 5.00E-07,
                                        'EIGWRITE': 0,
                                        }

        predictive_analysis = {'NPREDMAXMIN': 1,
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
        """Generate the PEST parameter file.

        For column 'PARNAME' the maximum length is 12 characters
        PARTRANS options = ['log', 'fixed', 'tied']
        PARCHGLIM options = ['factor', 'relative']

        Outputs a copy of the generated DataFrame as a csv if `csv_copy` property is set to `True`.

        :param params: dict, of PEST parameters. (Default value = None)

        :returns: dict, PEST parameters
        """

        header = ['PARNAME', 'PARTRANS', 'PARCHGLIM', 'PARVAL1', 'PARLBND', 'PARUBND',
                  'PARGP', 'SCALE', 'OFFSET', 'PARTIED', 'models', 'unit', 'comment']

        param_vals = list(params.values())
        num_param = len(list(params.keys()))
        PESTpar = pd.DataFrame(columns=header, index=list(params.keys()))
        PESTpar['PARNAME'] = list(params.keys())
        PESTpar['PARTRANS'] = [x.get('PARTRANS', 'log') for x in param_vals]
        PESTpar['PARCHGLIM'] = [x.get('PARCHGLIM', 'factor') for x in param_vals]
        PESTpar['PARVAL1'] = [x['PARVAL1'] for x in param_vals]
        PESTpar['PARLBND'] = [x.get('PARLBND', x['PARVAL1'] * 0.9) for x in param_vals]
        PESTpar['PARUBND'] = [x.get('PARUBND', x['PARVAL1'] * 1.1) for x in param_vals]
        PESTpar['PARGP'] = [x.get('PARGP', 'default') for x in param_vals]
        PESTpar['SCALE'] = [x.get('SCALE', 1.0) for x in param_vals]
        PESTpar['OFFSET'] = [x.get('OFFSET', 0.0) for x in param_vals]
        PESTpar['PARTIED'] = [''] * num_param
        PESTpar['models'] = [self.models_ID[0]] * num_param
        PESTpar['unit'] = ['-'] * num_param

        if self.csv_copy:
            self._createCSVcopyFromDataFrame(PESTpar, 'PESTpar')
        # end if

        return PESTpar
    # End PESTpar()

    def genPESTpgp(self):
        """Generate PEST parameter groups.

        Outputs a copy of the generated DataFrame as a csv if `csv_copy` property is set to `True`.

        :returns: None, `PESTgrp` key is updated in the `PEST_data` property"""

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
        PESTpgp['DERINC'] = [DERINCdefault] * num_parameter_groups
        PESTpgp['DERINCLB'] = [DERINCLBdefault] * num_parameter_groups
        PESTpgp['FORCEN'] = [FORCENdefault] * num_parameter_groups
        PESTpgp['DERINCMUL'] = [DERINCMULdefault] * num_parameter_groups
        PESTpgp['DERMTHD'] = [DERMTHDdefault] * num_parameter_groups

        if self.csv_copy:
            self._createCSVcopyFromDataFrame(PESTpgp, 'PESTpgp')
        # end if
        self.PEST_data['PESTpgp'] = PESTpgp
    # End genPESTpgp()

    def PESTobs(self, obs_grp):
        """Collate PEST observations into a DataFrame.

        Outputs a copy of the generated DataFrame as a csv if `csv_copy` property is set to `True`.

        :param obs_grp: dict, of observation groupings

        :returns: DataFrame, observation group time series
        """

        for index, key in enumerate(obs_grp.keys()):
            obs_grp_ts = obs_grp[key]['time_series'].copy()
            # Filter out null observations that were out of date range or
            # didn't map to the model mesh
            obs_grp_ts = obs_grp_ts[obs_grp_ts['obs_map'] != 'null']
            obs_grp_ts['OBGNME'] = obs_grp_ts['zone'] if obs_grp[key]['by_zone'] else key
            obs_grp_ts['WEIGHT'] = obs_grp[key]['weights']
            obs_grp_ts['model'] = self.models_ID[0]
            obs_grp_ts.rename(columns={'obs_map': 'OBSNME', 'value': 'OBSVAL'}, inplace=True)
            obs_grp_ts = obs_grp_ts[['OBGNME', 'OBSVAL', 'OBSNME', 'WEIGHT', 'model']]

            if index == 0:
                PESTobs = obs_grp_ts
            else:
                PESTobs = PESTobs.append(obs_grp_ts)  # DataFrame appends and returns new object
            # End if
        # End for

        if self.csv_copy:
            self._createCSVcopyFromDataFrame(PESTobs, 'PESTobs')
        # End if

        return PESTobs
    # End PESTobs()

    def _createCSVcopyFromDataFrame(self, df, name):
        """
        :param df:
        :param name:
        """

        fname = os.path.join(self.directory, name + '.csv')
        if os.path.exists(fname):
            print((fname + ' file exists already'))
        else:
            df.to_csv(fname, index=False)
        # end if
    # End _createCSVcopyFromDataFrame()

    def writePESTdict2csv(self):
        """TODO: Docs"""

        with open(os.path.join(self.directory, self.name + '.csv'), 'w') as f:
            w = csv.DictWriter(f, list(self.PEST_data.keys()))
            w.writeheader()
            w.writerow(self.PEST_data)
        # End with
    # End writePESTdict2csv()

    def genParameters(self, method='dataframe'):
        """Generate *name*_parameters.txt containing all models parameters.

        Parameters names and values are taken from 'PARNAME' and 'PARVAL1' in the PESTpar function.

        If method is 'csv', read in PEST parameters and observations from `PESTpar.csv` and `PESTobs.csv`

        :param method: str, one of 'csv', 'excel', or 'dataframe'. (Default value = 'dataframe')
        """

        # Ensure valid method argument is used
        assert method.lower() in ['csv', 'excel', 'dataframe'], "Method '{}' not recognized".format(method)

        print(('Generating {}_parameters.txt using {}'.format(self.name, method)))
        if (method == 'csv'):
            self.PEST_data['PESTpar'] = pd.read_csv(os.path.join(self.directory, 'PESTpar.csv'))
            self.PEST_data['PESTobs'] = pd.read_csv(os.path.join(self.directory, 'PESTobs.csv'))
        # End if

        # Generate *name*_parameters.txt
        self.PEST_data['PESTpar'].to_csv(os.path.join(self.directory + 'parameters.txt'), sep='\t',
                                         columns=['PARNAME', 'PARVAL1'], index=False)

    def genPestfiles(self, models_ID=None):
        """Generate a PEST folder and PEST files for a given chain of models_ID.

        All necessary PEST inputs are taken from PEST dict containing dicts or dataframes:
        * control data are read from PESTcon
        * parameters to calibrate are read from PESTpar, including PARTIED information
        * parameter groups are read from PESTpgp
        * observations are read from PESTobs
        * observation groups are taken as the unique entries of column OBGNME in PESTobs

        :param models_ID: list, of model IDs. Defaults to ['default']. (Default value = None)
        """

        print(('# Generating PEST files, %s #\n' % (datetime.datetime.now())))
        models_ID = ['default'] if not models_ID else models_ID
        PEST_name = 'pest'
        PEST_folder_name = '_'.join(['PEST'] + models_ID)

        PEST_folder = self.directory + os.path.sep
        if not os.path.isdir(PEST_folder):
            os.mkdir(PEST_folder)
        # end if

        # Retain only parameters corresponding to current models (uses the column 'models' of Dataframe PESTpar) by assigning 'fixed' to the other ones
        for row in self.PEST_data['PESTpar'].iterrows():
            if row[1]['models'] not in models_ID:
                self.PEST_data['PESTpar'].set_value(row[0], 'PARTRANS', 'fixed')
        # End for

        # Retain only observations corresponding to current models (uses the column 'model' of Dataframe PESTobs)
        required = []
        for row in self.PEST_data['PESTobs'].iterrows():
            if row[1]['model'] in models_ID:
                required += [row[0]]
        # end for

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
                print(('  ** prediction mode was detected\n' +
                      '  --> an observation ''prediction'' was automatically added that belongs to a group called ''predict'', conformally with what PEST requires\n' +
                      '  --> make sure to write the computed prediction at the last line of the observation file which must be writen after your model run\n'))
            # end if
        # end if

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

        PRIOR = {}

        # Generate PEST control file

        # Command line that PEST excutes to run the model
        # This assumes run.bat exists locally and contains a command like 'python run_and_postprocess_model.py'
        # which further assumes that run_and_postprocess_model.py exists and:
        # - loads the model and updates the parameters based on parameters.txt
        # - runs the model
        # - post-processes model results and writes relevant model outputs for PEST to read

        if os.name == 'nt':
            PESTCMD = 'run.bat'
        elif os.name == 'posix':
            PESTCMD = './run.sh'

        # Model input file
        INFLE = 'parameters.txt'

        # Corresponding template file for PEST to know how to write it
        TEMPFLE = 'parameters.tpl'

        OUTFLE = {}
        INSFLE = {}

        os.chdir(self.directory)
        for obs_gp in self.PEST_data['PESTobs']['OBGNME'].unique():
            # Model observation file (must be created by post-processing of model outputs)
            OUTFLE[obs_gp] = '.' + os.path.sep + 'model_' + models_ID[0] + os.path.sep + \
                'observations_' + obs_gp + '.txt'

            # Corresponding instruction file for PEST to know how to read it
            INSFLE[obs_gp] = 'observations_' + obs_gp + '.ins'
        # end for

        # PEST control file
        PESTFILE = PEST_name + '.pst'

        # Counters
        NPAR = self.PEST_data['PESTpar'].count()['PARNAME']
        NOBS = self.PEST_data['PESTobs'].count()['OBSNME']
        NPARGP = self.PEST_data['PESTpgp'].count()['PARGPNME']
        NPRIOR = len(PRIOR)
        NOBSGP = len(PESTobsgp)
        NTPLFILE = 1
        NINSFLE = len(INSFLE)

        # Open file
        with open(os.path.join(PEST_folder, PESTFILE), 'w') as f:
            c_data = self.PEST_data['PESTcon']['control_data']

            # SVD
            svd = self.PEST_data['PESTcon']['singular_value_decomposition']
            svd['MAXSING'] = len(self.PEST_data['PESTpar'].index)

            pcf_lines = [
                ['pcf'],
                ['* control data'],
                [c_data['RSTFLE'], c_data['PESTMODE']],
                [NPAR, NOBS, NPARGP, NPRIOR, NOBSGP],
                [NTPLFILE, NINSFLE, c_data['PRECIS'], c_data['DPOINT']],
                [c_data['RLAMBDA1'], c_data['RLAMFAC'], c_data['PHIRATSUF'], c_data['PHIREDLAM'], c_data['NUMLAM'],
                 c_data['JACUPDATE'], c_data['LAMFORGIVE'], c_data['DERFORGIVE']],
                [c_data['RELPARMAX'], c_data['FACPARMAX'], c_data['FACORIG']],
                [c_data['PHIREDSHW'], c_data['NOPTSWITCH'], c_data['BOUNDSCALE']],
                [c_data['NOPTMAX'], c_data['PHIREDSTP'], c_data['NPHISTP'], c_data['NPHINORED'], c_data['RELPARSTP'],
                 c_data['NRELPAR']],
                [c_data['ICOV'], c_data['ICOR'], c_data['IEIG'], c_data['IRES'], c_data['JCOSAVE'],
                 c_data['VERBOSEREC'], c_data['JCOSAVEITN'], c_data['REISAVEITN'], c_data['PARSAVEITN']]
            ]
            write_multiline(f, pcf_lines)

            svd_lines = [
                ['* singular value decomposition'],  # SVD
                [svd['SVDMODE']],
                [svd['MAXSING'], svd['EIGTHRESH']],
                [svd['EIGWRITE']]
            ]
            write_multiline(f, svd_lines)

            # parameter groups
            write_line(f, '* parameter groups')
            for row in self.PEST_data['PESTpgp'].iterrows():
                write_line(f, [str(x) for x in row[1].tolist()], delimit='\t')
            # end for

            # Parameter data
            write_line(f, '* parameter data')
            target_names = ['PARNAME', 'PARTRANS', 'PARCHGLIM', 'PARVAL1', 'PARLBND', 'PARUBND',
                            'PARGP', 'SCALE', 'OFFSET']
            for row in self.PEST_data['PESTpar'][target_names].iterrows():
                write_line(f, [str(x) for x in row[1].tolist()], delimit='\t')
            # end for
            for row in self.PEST_data['PESTpar'][['PARNAME', 'PARTRANS', 'PARTIED']].iterrows():
                if row[1]['PARTRANS'] == 'tied':
                    write_line(f, [row[1]['PARNAME'], row[1]['PARTIED']], delimit='\t')
                # end if
            # end for

            # Observation groups
            write_line(f, '* observation groups')
            write_multiline(f, [[str(x)] for x in PESTobsgp])

            # Observation data
            write_line(f, '* observation data')
            for row in self.PEST_data['PESTobs'][['OBSNME', 'OBSVAL', 'WEIGHT', 'OBGNME']].iterrows():
                write_line(f, [str(x) for x in row[1]], delimit='\t')
            # end for

            # Command line that pest executes
            write_line(f, '* model command line')
            write_line(f, [PESTCMD])

            # Model input/output
            write_line(f, '* model input/output')
            write_line(f, [TEMPFLE, INFLE], delimit='\t')
            for obs_gp in self.PEST_data['PESTobs']['OBGNME'].unique():
                write_line(f, [INSFLE[obs_gp], OUTFLE[obs_gp]], delimit='\t')
            # end for

            # PRIOR rules
            write_line(f, '* prior information')

            # Predictive analysis
            if self.predictive_analysis:
                pred_an = self.PEST_data['PESTcon']['predictive_analysis']
                prior_lines = [
                    ['* predictive analysis'],
                    [pred_an['NPREDMAXMIN']],
                    [pred_an['PD0'], pred_an['PD1'], pred_an['PD2']],
                    [pred_an['ABSPREDLAM'], pred_an['RELPREDLAM'], pred_an['INITSCHFAC'], pred_an['MULSCHFAC'],
                     pred_an['NSEARCH']],
                    [pred_an['ABSPREDSWH'], pred_an['RELPREDSWH']],
                    [pred_an['NPREDNORED'], pred_an['ABSPREDSTP'], pred_an['RELPREDSTP'], pred_an['NPREDSTP']]
                ]
                write_multiline(f, prior_lines)
        # end with

        # Generate initial parameters file
        # Generate PEST template file
        with open(os.path.join(PEST_folder, TEMPFLE), 'w') as f:
            write_line(f, 'ptf #')
            write_line(f, ['PARNAME', 'PARVAL'], delimit='\t')
            for row in self.PEST_data['PESTpar'][['PARNAME']].iterrows():
                f.write('%s\t#%-15s#\n' % (row[1]['PARNAME'], row[1]['PARNAME']))
                # write_line(f, [row[1]['PARNAME'], '{:15}'.format(row[1]['PARNAME'])], delimit='\t')
            # end for
        # end with

        # Generate PEST instruction file
        for obs_gp in self.PEST_data['PESTobs']['OBGNME'].unique():
            obs_pd = self.PEST_data['PESTobs']
            PESTobs_filtered = obs_pd[obs_pd['OBGNME'] == obs_gp]

            with open(os.path.join(PEST_folder, INSFLE[obs_gp]), 'w') as f:
                write_line(f, 'pif %%')
                for row in PESTobs_filtered.iterrows():
                    f.write('l1 !%s!\n' % (row[1]['OBSNME']))
                # end for
            # end with
        # end for

        # Generate a parameter uncertainty file using bound values
        # Calculate parameter standard deviation assuming normal distribution and that lower and upper bounds are 95% intervals
        STD = []
        for row in self.PEST_data['PESTpar'].iterrows():
            # If we assume that the upper and lower bounds specified represent
            # the lower and upper quartiles (i.e. the 25th and 75th percentiles)
            # of the likely parameter distribution,
            # which are + and - 2*sigma, then we can calculate the standard
            # deviation as 1/4 ( upper - lower) = 1/4 ( mean + 2sigma - (mean - 2sigma))
            # = 1/4 (4sigma) = sigma
            row_one = row[1]
            if row_one['PARTRANS'] == 'log':
                STD += [0.25 * (np.log10(row_one['PARUBND']) - np.log10(row_one['PARLBND']))]
            else:
                STD += [0.25 * (row_one['PARUBND'] - row_one['PARLBND'])]
            # end if
        # end for
        self.PEST_data['PESTpar']['STD'] = STD

        UNCERTAINTYFILE = PEST_name + '.unc'
        with open(os.path.join(PEST_folder, UNCERTAINTYFILE), 'w') as f:
            write_line(f, '# Parameter uncertainty file')
            write_line(f, '# for filling C(k)\n')  # new line is intentional

            # Needs to be for every other parameter
            write_line(f, 'START STANDARD_DEVIATION')
            for row in self.PEST_data['PESTpar'].iterrows():
                if row[1]['PARTRANS'] != 'fixed':
                    write_line(f, [row[1]['PARNAME'], row[1]['STD']], delimit='\t')
                # end if
            # end for
            write_line(f, 'END STANDARD_DEVIATION')
        # end with

        print('\nPEST files generated\n')
        print((' %s\n' % PESTFILE))
        print('\nPEST files generation completed!\n')

    def updateparameterswithpestbestpar(self, pestparfile):
        """Generate parameters.txt containing all models parameters.

        Parameters names and values are taken from 'PARNAME' and 'PARVAL1' in PESTpar

        :param pestparfile:
        """

        raise NotImplementedError("`updateparameterswithpestbestpar` method is incomplete and should not be used.")
        print('Generating parameters.txt\n')

        # Check existence of the PEST .par file and read it
        if not os.path.exists(pestparfile):
            raise IOError("Can't find PEST .par file <<{}>>".format(pestparfile))
        # end if

        # with open(pestparfile, 'r') as f:
        #     text = f.readlines()
        # # end with

        # PESTbestpar = textscan(fileID,'%s %f %f %f\n','HeaderLines',1,'Delimiter',' ','MultipleDelimsAsOne',1);
        # if isnan(PESTbestpar{1,2})
        #     fclose(fileID);
        #     fileID = fopen(pestparfile);
        #     PESTbestpar = textscan(fileID,'%s %f %f %f','HeaderLines',1,'Delimiter','\t','MultipleDelimsAsOne',1);
        # # end

        # Assign PEST best parameter value

        self.PEST_data['PESTpar']['PARVAL1'] = PESTbestpar

        # Generate *name*_parameters.txt

        self.PEST_data['PESTpar'].to_csv(self.directory + 'best_parameters.txt',
                                         sep='\t', columns=['PARNAME', 'PARVAL1'], index=False)
