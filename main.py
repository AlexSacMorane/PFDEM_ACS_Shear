# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the main file.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import os
import shutil
from pathlib import Path

#Own function and class
import Create_IC
import Create_IC_Polygonal
import Confine_Polygonal
import Contact_gg
import Contact_gimage
import Etai
import Grain
import Shear_Polygonal
import Report
import User
import Owntools
import Owntools.Plot
import Owntools.Save

#-------------------------------------------------------------------------------
#Functions
#-------------------------------------------------------------------------------

def plan_simulation():
    """
    Plan the simulation :
        - Get data
        - Create folders and files
        - Find a name for the simulation

        Input :
            Nothing
        Output :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a simulation report (a report)
    """
    #get data
    dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations = User.All_parameters()

    #debug folder
    if Path('Debug').exists():
        shutil.rmtree('Debug')
    os.mkdir('Debug')
    if dict_algorithm['Debug'] or dict_algorithm['Debug_DEM'] or dict_ic['Debug_DEM_IC'] :
        if dict_algorithm['Debug_DEM'] :
            os.mkdir('Debug/Shear')
        if dict_ic['Debug_DEM'] :
            os.mkdir('Debug/Init_disks')
            os.mkdir('Debug/Init_polygons')
            os.mkdir('Debug/Init_polygons_group')

    #report
    simulation_report = Report.Report('Debug/Report.txt')

    #find name
    if dict_algorithm['SaveData'] :
        #check if save folder exists
        if not Path('../'+dict_algorithm['main_folder_name']).exists():
            shutil.mkdir('../'+dict_algorithm['main_folder_name'])
        i_run = 1
        folderpath = Path('../'+dict_algorithm['main_folder_name']+'/'+dict_algorithm['template_simulation_name']+str(i_run))
        while folderpath.exists():
            i_run = i_run + 1
            folderpath = Path('../'+dict_algorithm['main_folder_name']+'/'+dict_algorithm['template_simulation_name']+str(i_run))
        dict_algorithm['name_folder'] = dict_algorithm['template_simulation_name']+str(i_run)
    else :
        dict_algorithm['name_folder'] = 'One_Run'

    return dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report

#-------------------------------------------------------------------------------

def generate_ic(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report) :
    """
    Generate an initial condition.

    First, disks are generated and loaded. Once a steady-state is detected, grains are discretized.
    Then, polygons are loaded until a steady-state is detected.

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a simulation report (a report)
        Output :
            Nothing, but dictionnaries are updated
    """
    #generate and load disks
    simulation_report.tic_tempo()
    Create_IC.LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)
    simulation_report.tac_tempo('Loading with disks')

    #discretization
    simulation_report.write_and_print('Discretize the sample\n', 'Discretize the sample\n')
    simulation_report.tic_tempo()
    dict_ic = Create_IC_Polygonal.Discretize_Grains(dict_ic, dict_geometry['discretization']) #overwrite sphere dict_ic
    simulation_report.tac_tempo('From disks to polygons')

    #load discrete grains
    print('Load with top plate')
    simulation_report.tic_tempo()
    Create_IC_Polygonal.DEM_loading(dict_algorithm, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)
    simulation_report.tac_tempo('Loading with polygons')

#-------------------------------------------------------------------------------

def define_group(dict_algorithm, dict_ic, dict_sample, simulation_report):
    """
    Define top and bottom groups and walls are deleted.

        Input :
             an algorithm dictionnary (a dict)
             an initial dictionnary (a dict)
             a sample dictionnary (a dict)
             a simulation report (a report)
         Output :
            Nothing, but dictionnaries are updated
    """
    #define packs
    simulation_report.tic_tempo()
    i_bottom = 0
    i_top = 0
    for grain in dict_ic['L_g_tempo'] :
        if grain.is_group(dict_sample['y_box_min'], dict_sample['y_box_min'] + dict_algorithm['bottom_height'], 'Bottom') :
            i_bottom = i_bottom + 1
        elif grain.is_group(dict_sample['y_box_max'] - dict_algorithm['top_height'], dict_sample['y_box_max'], 'Top') :
            i_top = i_top + 1
    simulation_report.write_and_print('Define groups\n','Define groups')
    simulation_report.write_and_print(str(i_bottom)+' grains in Bottom group\n'+str(i_top)+' grains in Top group\n\n', str(i_bottom)+' grains in Bottom group\n'+str(i_top)+' grains in Top group\n')

    #delete contact gw
    dict_ic['L_contact_gw'] = []
    dict_ic['L_contact_gw_ij'] = []

    #plot group distribution
    if dict_algorithm['Debug']:
        Owntools.Plot.Plot_group_distribution(dict_ic)
        simulation_report.tac_tempo('Define groups')

#-------------------------------------------------------------------------------

def load_ic_group(dict_algorithm, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report) :
    """
    Generate an initial condition.

    Polygons are loaded with top group until a steady-state is detected.

        Input :
            an algorithm dictionnary (a dict)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a simulation report (a report)
        Output :
            Nothing, but dictionnaries are updated
    """
    #load discrete grains
    print('Load with top group')
    simulation_report.tic_tempo()
    Create_IC_Polygonal.DEM_loading_group(dict_algorithm, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)
    simulation_report.tac_tempo('Loading with polygons with groups')

#-------------------------------------------------------------------------------

def from_ic_to_real(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report):
    """
    Convert initial configuration to a current one.

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a simulation report (a report)
        Output :
            Nothing, but dictionnaries are updated
    """
    simulation_report.write_and_print('Convert initial configuration to current one\n\n', 'Convert initial configuration to current one\n')
    Owntools.convert_ic_to_real(dict_ic, dict_sample)
    #save
    Owntools.Save.save_dicts_ic(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)

#-------------------------------------------------------------------------------

def load_sample(dict_algorithm, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report):
    '''
    Load the sample.

        Input :
            an algorithm dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a tracker dictionnary (a dict)
            a simulation report (a report)
        Output :
            Nothing, but dictionnaries are updated
    '''
    #shear sample
    simulation_report.write_and_print('\nConfine the sample\n', 'Confine the sample')
    simulation_report.tic_tempo()
    Confine_Polygonal.DEM_confine_load(dict_algorithm, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report)
    simulation_report.tac_tempo('Confinement')

#-------------------------------------------------------------------------------

def define_etai(dict_algorithm, dict_material, dict_sample, simulation_report):
    '''
    Define the etai to describe the sample.

        Input :
            an algorithm dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a simulation report (a report)
        Output :
            Nothing, but sample dictionnary is updated with etais
    '''
    simulation_report.write_and_print('\nDefining the etais\n', 'Defining the etais')
    simulation_report.tic_tempo()
    #build phase field for grains
    for grain in dict_sample['L_g']:
        grain.build_etai_M(dict_material,dict_sample)
    #build etais
    dict_sample['L_etai'] = list([])
    for group in ['Top', 'Current', 'Bottom'] :
        Etai.etai_distribution(dict_algorithm, dict_sample, simulation_report, group)
    if dict_algorithm['Debug']:
        Owntools.Plot.Plot_etais(dict_sample)
    simulation_report.write_and_print(f"{round(len(dict_sample['L_g'])/len(dict_sample['L_etai']),1)} grains in average per etai.\n\n", f"{round(len(dict_sample['L_g'])/len(dict_sample['L_etai']),1)} grains in average per etai.\n")
    simulation_report.tac_tempo('Defining the etais')

#-------------------------------------------------------------------------------

def shear_sample(dict_algorithm, dict_geometry, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report):
    '''
    Shear the sample.

        Input :
            an algorithm dictionnary (a dict)
            a geometry disctionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a tracker dictionnary (a dict)
            a simulation report (a report)
        Output :
            Nothing, but dictionnaries are updated
    '''
    #shear sample
    simulation_report.write_and_print('\nShearing the sample\n', 'Shearing the sample')
    simulation_report.tic_tempo()
    Shear_Polygonal.DEM_shear_load(dict_algorithm, dict_geometry, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report)
    simulation_report.write_and_print('\n','')
    simulation_report.tac_tempo('Shearing')

#-------------------------------------------------------------------------------

def close_simulation(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report):
    '''
    Close the simulation.

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a tracker dictionnary (a dict)
            a simulation report (a report)
        Output :
            Nothing, but dictionnaries are updated
    '''
    #make movie of the different configuration
    if dict_algorithm['Debug_DEM'] :
        Owntools.Plot.Plot_mp4('Debug/Shear/Config_','Debug/Shear/Configuration.mp4')

    simulation_report.end()

    #if dict_algorithm['cleanData'] :
    #    pass

    #final save
    if dict_algorithm['SaveData']:
        Owntools.Save.save_dicts(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report)
        name_actual_folder = os.path.dirname(os.path.realpath(__file__))
        shutil.copytree(name_actual_folder, '../'+dict_algorithm['main_folder_name']+'/'+dict_algorithm['name_folder'])
        os.remove(dict_algorithm['name_folder']+'_save_dicts')
        os.remove(dict_algorithm['name_folder']+'_ic_save_dicts')

#-------------------------------------------------------------------------------
#main
#-------------------------------------------------------------------------------

if '__main__' == __name__:

    #open simulation and get data
    dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report = plan_simulation()

    #ic generation
    generate_ic(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)
    define_group(dict_algorithm, dict_ic, dict_sample, simulation_report)

    #convert ic to real grain
    from_ic_to_real(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)

    #create mesh
    Owntools.create_mesh(dict_algorithm, dict_material, dict_sample)

    #create etais
    define_etai(dict_algorithm, dict_material, dict_sample, simulation_report)

    #add solute
    User.generate_solute(dict_sample)

    #load
    dict_tracker = {}
    shear_sample(dict_algorithm, dict_geometry, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report)

    #close simulation
    close_simulation(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report)
