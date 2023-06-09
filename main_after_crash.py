# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the file to restart a simulation after a crash.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import os
import pickle
import shutil
from pathlib import Path

#Own function and class
import main
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
#User
#-------------------------------------------------------------------------------

name_to_load = 'Run_9_ic_save_dicts'

#-------------------------------------------------------------------------------
#load data
#-------------------------------------------------------------------------------

toload = open(name_to_load,'rb')
dict_save = pickle.load(toload)
toload.close()
dict_algorithm = dict_save['algorithm']
dict_geometry = dict_save['geometry']
dict_ic = dict_save['ic']
dict_material = dict_save['material']
dict_sample = dict_save['sample']
dict_sollicitations = dict_save['sollicitations']
simulation_report = dict_save['report']

#-------------------------------------------------------------------------------
#plan simulation
#-------------------------------------------------------------------------------

simulation_report.write('\nA crash occurs...\n\n')

#-------------------------------------------------------------------------------
#main
#-------------------------------------------------------------------------------

if name_to_load[-14:] =='_ic_save_dicts':
    #create mesh
    Owntools.create_mesh(dict_algorithm, dict_material, dict_sample)

    #create etais
    main.define_etai(dict_algorithm, dict_material, dict_sample, simulation_report)

    #add solute
    User.generate_solute(dict_sample)

    #load
    dict_tracker = {}
    main.shear_sample(dict_algorithm, dict_geometry, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report)

    #close simulation
    main.close_simulation(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report)
