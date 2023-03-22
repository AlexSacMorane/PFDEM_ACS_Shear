# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions to save data.
"""

#-------------------------------------------------------------------------------
#Librairies
#-------------------------------------------------------------------------------

import pickle

#Own
import Contact_gg
import Contact_gimage
import Grain
import Report

#-------------------------------------------------------------------------------
#functions
#-------------------------------------------------------------------------------

def save_dicts_debug(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report):
    """
    Save dictionnaries after the initial configuration.

    To be deleted. Only for debug

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a simulation report (a report)
        Output :
            Nothing, but a save file is generated (a file)
    """
    outfile = open(dict_algorithm['name_folder']+'_db_save_dicts','wb')
    dict_save = {}
    dict_save['algorithm'] = dict_algorithm
    dict_save['geometry'] = dict_geometry
    dict_save['ic'] = dict_ic
    dict_save['material'] = dict_material
    dict_save['sample'] = dict_sample
    dict_save['sollicitations'] = dict_sollicitations
    dict_save['report'] = simulation_report
    pickle.dump(dict_save,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_dicts_ic(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report):
    """
    Save dictionnaries after the initial configuration.

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a simulation report (a report)
        Output :
            Nothing, but a save file is generated (a file)
    """
    outfile = open(dict_algorithm['name_folder']+'_ic_save_dicts','wb')
    dict_save = {}
    dict_save['algorithm'] = dict_algorithm
    dict_save['geometry'] = dict_geometry
    dict_save['ic'] = dict_ic
    dict_save['material'] = dict_material
    dict_save['sample'] = dict_sample
    dict_save['sollicitations'] = dict_sollicitations
    dict_save['report'] = simulation_report
    pickle.dump(dict_save,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_dicts(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report):
    """
    Save dictionnaries at the end of PFDEM iteration.

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a save file is generated (a file)
    """
    outfile = open(dict_algorithm['name_folder']+'_save_dicts','wb')
    dict_save = {}
    dict_save['algorithm'] = dict_algorithm
    dict_save['geometry'] = dict_geometry
    dict_save['ic'] = dict_ic
    dict_save['material'] = dict_material
    dict_save['sample'] = dict_sample
    dict_save['sollicitations'] = dict_sollicitations
    dict_save['tracker'] = dict_tracker
    dict_save['report'] = simulation_report
    pickle.dump(dict_save,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_dicts_before_pf(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report):
    """
    Save dictionnaries at the end of PFDEM iteration.

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a save file is generated (a file)
    """
    outfile = open(dict_algorithm['name_folder']+'_save_dicts_before_pf','wb')
    dict_save = {}
    dict_save['algorithm'] = dict_algorithm
    dict_save['geometry'] = dict_geometry
    dict_save['ic'] = dict_ic
    dict_save['material'] = dict_material
    dict_save['sample'] = dict_sample
    dict_save['sollicitations'] = dict_sollicitations
    dict_save['tracker'] = dict_tracker
    dict_save['report'] = simulation_report
    pickle.dump(dict_save,outfile)
    outfile.close()
