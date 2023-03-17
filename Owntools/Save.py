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
import Contact
import Contact_gw
import Grain
import Report

#-------------------------------------------------------------------------------
#functions
#-------------------------------------------------------------------------------

def save_dicts(dict_algorithm, dict_geometry, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report):
    """
    Save dictionnaries at the end of PFDEM iteration.

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
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
    dict_save['material'] = dict_material
    dict_save['sample'] = dict_sample
    dict_save['sollicitations'] = dict_sollicitations
    dict_save['tracker'] = dict_tracker
    dict_save['report'] = simulation_report
    pickle.dump(dict_save,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_dicts_before_pf(dict_algorithm, dict_geometry, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report):
    """
    Save dictionnaries at the end of PFDEM iteration.

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
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
    dict_save['material'] = dict_material
    dict_save['sample'] = dict_sample
    dict_save['sollicitations'] = dict_sollicitations
    dict_save['tracker'] = dict_tracker
    dict_save['report'] = simulation_report
    pickle.dump(dict_save,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_tempo(dict_algorithm,dict_tracker):
    """
    Save trackers and configuration during  PFDEM iteration.

        Input :
            an algorithm dictionnary (a dict)
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a save file is generated (a file)
    """
    outfile = open('../'+dict_algorithm['main_folder_name']+'/'+dict_algorithm['name_folder']+'_save_tempo','wb')
    dict_save = {}
    dict_save['k0_xmin_L'] = dict_tracker['k0_xmin_L']
    dict_save['k0_xmax_L'] = dict_tracker['k0_xmax_L']
    dict_save['S_dissolved_L'] = dict_tracker['S_dissolved_L']
    dict_save['S_grains_dissolvable_L'] = dict_tracker['S_grains_dissolvable_L']
    dict_save['S_dissolved_perc_L'] = dict_tracker['S_dissolved_perc_L']
    dict_save['S_dissolved_perc_dissolvable_L'] = dict_tracker['S_dissolved_perc_dissolvable_L']
    dict_save['S_grains_L'] = dict_tracker['S_grains_L']
    pickle.dump(dict_save,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_final(dict_algorithm,dict_tracker):
    """
    Save trackers and configuration at the end of simulation.

        Input :
            an algorithm dictionnary (a dict)
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a save file is generated (a file)
    """
    os.remove('../'+dict_algorithm['main_folder_name']+'/'+dict_algorithm['name_folder']+'_save_tempo')
    outfile = open('../'+dict_algorithm['main_folder_name']+'/'+dict_algorithm['name_folder']+'_save','wb')
    dict_save = {}
    dict_save['k0_xmin_L'] = dict_tracker['k0_xmin_L']
    dict_save['k0_xmax_L'] = dict_tracker['k0_xmax_L']
    dict_save['S_dissolved_L'] = dict_tracker['S_dissolved_L']
    dict_save['S_dissolved_perc_L'] = dict_tracker['S_dissolved_perc_L']
    dict_save['S_grains_L'] = dict_tracker['S_grains_L']
    pickle.dump(dict_save,outfile)
    outfile.close()
