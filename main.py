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
from datetime import datetime
from pathlib import Path

#Own function and class
import Create_IC
import Create_IC_Polygonal
import Confine_Polygonal
import Contact_gg
import Contact_gimage
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
    simulation_report = Report.Report('Debug/Report.txt',datetime.now())

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
    simulation_report.tic_tempo(datetime.now())
    Create_IC.LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)
    simulation_report.tac_tempo(datetime.now(), 'Loading with disks')

    #discretization
    simulation_report.write_and_print('Discretize the sample\n', 'Discretize the sample\n')
    simulation_report.tic_tempo(datetime.now())
    dict_ic = Create_IC_Polygonal.Discretize_Grains(dict_ic, dict_geometry['discretization']) #overwrite sphere dict_ic
    simulation_report.tac_tempo(datetime.now(), 'From disks to polygons')

    #load discrete grains
    print('Load with top plate')
    simulation_report.tic_tempo(datetime.now())
    Create_IC_Polygonal.DEM_loading(dict_algorithm, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)
    simulation_report.tac_tempo(datetime.now(), 'Loading with polygons')

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
    simulation_report.tic_tempo(datetime.now())
    i_bottom = 0
    i_top = 0
    for grain in dict_ic['L_g_tempo'] :
        if grain.is_group(dict_sample['y_box_min'], dict_sample['y_box_min'] + dict_algorithm['bottom_height'], 'Bottom') :
            i_bottom = i_bottom + 1
        elif grain.is_group(dict_sample['y_box_max'] - dict_algorithm['top_height'], dict_sample['y_box_max'], 'Top') :
            i_top = i_top + 1
    simulation_report.write_and_print('\nDefine groups\n','\nDefine groups')
    simulation_report.write_and_print(str(i_bottom)+' grains in Bottom group\n'+str(i_top)+' grains in Top group\n\n', str(i_bottom)+' grains in Bottom group\n'+str(i_top)+' grains in Top group\n')

    #delete contact gw
    dict_ic['L_contact_gw'] = []
    dict_ic['L_contact_gw_ij'] = []

    #plot group distribution
    Owntools.Plot.Plot_group_distribution(dict_ic)
    simulation_report.tac_tempo(datetime.now(), 'Define groups')

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
    simulation_report.tic_tempo(datetime.now())
    Create_IC_Polygonal.DEM_loading_group(dict_algorithm, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)
    simulation_report.tac_tempo(datetime.now(), 'Loading with polygons with groups')

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
    simulation_report.write_and_print('\nConvert initial configuration to current one\n\n', '\nConvert initial configuration to current one\n')
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
    simulation_report.tic_tempo(datetime.now())
    Confine_Polygonal.DEM_confine_load(dict_algorithm, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report)
    simulation_report.tac_tempo(datetime.now(), 'Confinement')

#-------------------------------------------------------------------------------

def shear_sample(dict_algorithm, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report):
    '''
    Shear the sample.

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
    simulation_report.write_and_print('\nShearing the sample\n', 'Shearing the sample')
    simulation_report.tic_tempo(datetime.now())
    Shear_Polygonal.DEM_shear_load(dict_algorithm, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report)
    simulation_report.tac_tempo(datetime.now(), 'Shearing')

#-------------------------------------------------------------------------------

def close_simulation(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report):
    '''
    Close the simulation.

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
    '''
    pass

#-------------------------------------------------------------------------------
#main
#-------------------------------------------------------------------------------

if '__main__' == __name__:

    #open simulation and get data
    dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report = plan_simulation()

    #ic generation
    generate_ic(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)
    define_group(dict_algorithm, dict_ic, dict_sample, simulation_report)
    #load_ic_group(dict_algorithm, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)

    #convert ic to real grain
    from_ic_to_real(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)

    #load
    dict_tracker = {}

    load_sample(dict_algorithm, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report)
    raise ValueError('stop')

    shear_sample(dict_algorithm, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report)

    #cose simulation
    close_simulation(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, simulation_report)
