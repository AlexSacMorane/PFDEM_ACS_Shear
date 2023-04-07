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

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def pf_simulation():
    """
    Convert from DEM, prepare simulation, run simulation and convert to DEM.

        Input :

        Output :

    """
    #compute E_mec
    Owntools.Compute.Compute_Emec(dict_material, dict_sample, dict_sollicitations)

    #compute kc
    Owntools.Compute.Compute_kc(dict_algorithm, dict_material, dict_sample)

    #write etas

    #write solute

    #write kc

    #write emec

    #write .i

    #run simulation

    #sort files

    #pf to dem
