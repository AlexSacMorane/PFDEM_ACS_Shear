# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains functions used in the simulation to create and run phase field simulations.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import os
from pathlib import Path

#Own
import Owntools
import Owntools.Compute
import Owntools.Write

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def pf_simulation(dict_algorithm, dict_material, dict_sample, dict_sollicitations):
    """
    Convert from DEM, prepare simulation, run simulation and convert to DEM.

        Input :

        Output :

    """
    #plan simulation
    L_folder = ['Data', 'Input', 'Output']
    for folder in L_folder :
        if not Path(folder).exists():
            os.mkdir(folder)
        else :
            if dict_algorithm['clean_memory']:
                shutil.rmtree(folder)
                os.mkdir(folder)

    #compute data
    Owntools.Compute.Compute_Emec(dict_material, dict_sample, dict_sollicitations)
    Owntools.Compute.Compute_kc(dict_algorithm, dict_material, dict_sample)

    #write data
    Owntools.Write.Write_eta_txt(dict_algorithm, dict_sample)
    Owntools.Write.Write_solute_txt(dict_algorithm, dict_sample)
    Owntools.Write.Write_kc_txt(dict_algorithm, dict_sample)
    Owntools.Write.Write_Emec_txt(dict_algorithm, dict_sample)

    #write .i
    Owntools.Write.Write_i(dict_algorithm, dict_material, dict_sample, dict_sollicitations)

    #run simulation
    os.system('mpiexec -n '+str(dict_algorithm['np_proc'])+' ~/projects/moose/modules/combined/combined-opt -i '+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'.i')

    #sort files
    j_str = Owntools.Sort_Files(dict_algorithm)

    #pf to dem
