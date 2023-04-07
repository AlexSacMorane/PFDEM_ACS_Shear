# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains functions used in the simulation to create and run phase field simulations.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

#Own
import Owntools.Compute
import Owntools.Write

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def pf_simulation():
    """
    Convert from DEM, prepare simulation, run simulation and convert to DEM.

        Input :

        Output :

    """
    #compute data
    Owntools.Compute.Compute_Emec(dict_material, dict_sample, dict_sollicitations)
    Owntools.Compute.Compute_kc(dict_algorithm, dict_material, dict_sample)

    #write data
    Owntools.Write.Write_eta_txt(dict_algorithm, dict_sample)
    Owntools.Write.Write_solute_txt(dict_algorithm, dict_sample)
    Owntools.Write.Write_kc_txt(dict_algorithm, dict_sample)
    Owntools.Write.Write_Emec_txt(dict_algorithm, dict_sample)

    #write .i

    #run simulation

    #sort files

    #pf to dem
