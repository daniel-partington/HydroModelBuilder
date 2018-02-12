# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 10:15:07 2017

@author: part0075
"""

import subprocess


def PESTutil_run(name):
    """


    :param name: 

    """

    try:
        print("Running {}:\n".format(name))
        subprocess.check_output([name, '<', name + '.in'], shell=True)
    except subprocess.CalledProcessError as e:
        print("stdout output on error:\n" + e.output)


'''
FROM "pest_settings.pdf":

Here is a general “recipe” for PEST setup.

1. Build a PEST control, template and instruction files in the usual manner.
   In the PEST control file file:
     Set initial parameter values at their preferred values based on expert
      knowledge.
     Assign parameters of different types to different parameter groups. Make
      sure that the names of these groups are 6 characters or less in length.
      (This will allow the ADDREG1 utility to formulate unique names for
      regularization groups which constrain them – see below.)
     Assign appropriate bounds to parameters, these reflecting the endpoints
      of their range of realistic values.
     Assign observations of different types to different observation groups.
     Weight observations within each group in a way that provides suitable
      weighting relatively. For some groups this will mean that all observations
      are assigned a weight of 1.0. For other groups weights may need to reflect
      the values of the observations to which they pertain, especially where the
      latter can vary over many orders of magnitude, as can river flows and
      contaminant concentrations in groundwater. In such cases higher weights
      may need to be provided to low-valued observations so that they are visible to PEST.
     Log transform all parameters whose lower bound is greater than 0.0.
      (Provide a lower bound that is slightly greater than zero, rather than
      exactly zero, if necessary.)
     Set RLAMDA1 to 10.0 and RLAMFAC to -3; the latter guarantees fast movement
      of the Marquardt lambda.

2. Introduce a “singular value decomposition” section to the PEST control
   file. Use of singular value decomposition (SVD) ensures that PEST maintains
   numerical stability, regardless of how ill-posed is the inverse problem.
   (In its early days, PEST sometimes failed to work where inverse problems
   were ill-posed.) In the “singular value decomposition” section set:
     SVDMODE to 1;
     MAXSING to the number of adjustable parameters;
     EIGTHRESH to 5×10-7;
     EIGWRITE to 0.

3. However if you are estimating more than about 2500 parameters, introduce an
   “lsqr” section to the PEST control file instead of a “singular value
   decomposition” section. In this section set:
     LSQRMODE to 1;
     LSQR_ATOL and LSQR_BTOL to 10-6;
     LSQR_CONLIM to 1000;
     LSQR_ITNLIM to 4 times the number of adjustable parameters;
     LSQRWRITE to 0.

4. Set NOPTMAX to zero in the “control data” section of the PEST control file
   and then run PEST. PEST will run the model once and calculate an objective
   function, together with the contribution made to the total objective
   function by each observation group. This will allow you to ensure that all
   is working as it should. It will also provide information needed by the
   PWTADJ1 utility, which is run next.

5. Run PWTADJ1. This creates a new PEST control file in which weights are
   adjusted such that the contribution made to the total objective function by
   each observation group is the same. This prevents the information content of
   any group from being invisible to the inversion process.

6. Add Tikhonov regularization to the PEST control file by running the ADDREG1
   utility. ADDREG1 automatically sets the target measurement objective
   function (PHIMLIM) to 10-10. This, of course, is way too low; however it is
   often a good idea to attain the lowest measurement objective function that
   you can on the initial PEST run. Then, on later runs, set PHIMLIM to a value
   somewhat higher than the lowest achievable objective function in order to
   keep estimated parameter values closer to their initial values if you judge
   that this is a reasonable thing to do. Set PHIMACCEPT about 2% higher than
   PHIMLIM.

   In adding regularization to the PEST control file ADDREG1 automatically does
   the following.
     It provides a prior information equation for each parameter in which the
      preferred value of that parameter is equated to its initial value.
     It assigns these prior information equations to regularization groups
      that inherit their names from the parameter groups to which the
      parameters belong. (It must add “regul_” as a prefix to the name of each
      parameter group; that is why parameter group names should be 6 characters
      or less in length.)
     It sets the IREGADJ variable to 1. This allows PEST to assign different
      regularization weight factors to different prior information groups.
     It sets the FRACPHIM variable to 0.1. This instructs PEST to temporarily
      adopt a target measurement objective function of one tenth the current
      measurement objective function during all iterations, unless PHIMLIM is
      greater than this, in which case PHIMLIM prevails.

7. Run PEST.
'''

run_utilities = ['PWTADJ1', 'ADDREG1', ]

for util in run_utilities:
    PESTutil_run(util)
