# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains functions used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np
import random
import math

#Own
from Grain import Grain, Grain_Image
import Contact_gg
import Contact_gimage
import Owntools
import Owntools.Confinement
import Owntools.Plot

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def DEM_confine_load(dict_algorithm, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report):
    """
    Loading the granular system with vertical load.

        Input :
            an algorithm dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitation dictionnary (a dict)
            a tracker dictionnary (a dict)
            a simultion report (a report)
        Output :
            Nothing, but sample dictionnary is updated
    """
    dict_algorithm['i_DEM'] = 0
    i_DEM_0 = dict_algorithm['i_DEM']
    DEM_loop_statut = True
    #Initialisation
    dict_sample['L_contact'] = []
    dict_sample['L_contact_ij'] = []
    dict_sample['L_contact_gimage'] = []
    dict_sample['L_contact_ij_gimage'] = []
    dict_sample['id_contact'] = 0
    dict_sample['L_g_image'] = []
    dict_sample['L_i_image'] = []
    #trackers
    dict_tracker['vertical_force_L'] = []
    dict_tracker['compacity_L'] = []
    dict_tracker['dy_top_L'] = []
    dict_tracker['Ecin_L'] = []
    dict_tracker['Force_L'] = []

    while DEM_loop_statut :

        dict_algorithm['i_DEM'] = dict_algorithm['i_DEM'] + 1

        #create image
        for grain in dict_sample['L_g']:
            #left wall
            if (grain.center[0] - dict_sample['x_box_min']) < dict_algorithm['d_to_image'] :
                if grain.id in dict_sample['L_i_image'] : #image exists
                    image = dict_sample['L_g_image'][dict_sample['L_i_image'].index(grain.id)]
                    if image.position == 'right' :
                        image.position = 'left'
                else : #image does not exist
                    dict_sample['L_g_image'].append(Grain_Image(grain, 'left'))
                    dict_sample['L_i_image'].append(grain.id)
            #right wall
            elif (dict_sample['x_box_max'] - grain.center[0]) < dict_algorithm['d_to_image'] :
                if grain.id in dict_sample['L_i_image'] : #image exists
                    image = dict_sample['L_g_image'][dict_sample['L_i_image'].index(grain.id)]
                    if image.position == 'left' :
                        image.position = 'right'
                else : #image does not exist
                    dict_sample['L_g_image'].append(Grain_Image(grain, 'right'))
                    dict_sample['L_i_image'].append(grain.id)
            #center
            else :
                if grain.id in dict_sample['L_i_image'] : #image exists
                    i_toremove = dict_sample['L_i_image'].index(grain.id)
                    dict_sample['L_g_image'].pop(i_toremove)
                    dict_sample['L_i_image'].pop(i_toremove)
                    L_i_toremove = []
                    for ij_gimage in dict_sample['L_contact_ij_gimage'] :
                        if grain.id == ij_gimage[1] :
                            L_i_toremove.append(dict_sample['L_contact_ij_gimage'].index(ij_gimage))
                    L_i_toremove.reverse()
                    for i_toremove in L_i_toremove:
                        dict_sample['L_contact_gimage'].pop(i_toremove)
                        dict_sample['L_contact_ij_gimage'].pop(i_toremove)
        #translate image
        for image in dict_sample['L_g_image']:
            if image.position == 'left' :
                image.translation(np.array([dict_sample['x_box_max'] - dict_sample['x_box_min'], 0]))
            elif image.position == 'right' :
                image.translation(np.array([dict_sample['x_box_min'] - dict_sample['x_box_max'], 0]))

        #Contact detection
        if (dict_algorithm['i_DEM']-i_DEM_0-1) % dict_algorithm['i_update_neighborhoods']  == 0:
            Contact_gg.Update_Neighborhoods(dict_algorithm, dict_sample)
            Contact_gimage.Update_Neighborhoods(dict_algorithm, dict_sample)
        Contact_gg.Grains_contact_Neighborhoods(dict_sample,dict_material)
        Contact_gimage.Grains_contact_Neighborhoods(dict_sample,dict_material)

        #Sollicitation computation
        for grain in dict_sample['L_g']:
             grain.init_F_control(dict_sollicitations['gravity'])
        for contact in  dict_sample['L_contact']+dict_sample['L_contact_gimage']:
            #do not consider the contact inside top and bottom groups
            if not (contact.g1.group == 'Top' and contact.g2.group =='Top') or not (contact.g1.group == 'Bottom' and contact.g2.group =='Bottom') :
                contact.normal()
                contact.tangential(dict_algorithm['dt_DEM'])

        #Delete contacts gg and gimage with no overlap
        L_i_toremove = []
        for i_contact in range(len(dict_sample['L_contact'])):
            if dict_sample['L_contact'][i_contact].overlap_normal < 0:
                L_i_toremove.append(i_contact)
        L_i_toremove.reverse()
        for i_toremove in L_i_toremove:
            dict_sample['L_contact'].pop(i_toremove)
            dict_sample['L_contact_ij'].pop(i_toremove)
        L_i_toremove = []
        for i_contact in range(len(dict_sample['L_contact_gimage'])):
            if dict_sample['L_contact_gimage'][i_contact].overlap_normal < 0:
                L_i_toremove.append(i_contact)
        L_i_toremove.reverse()
        for i_toremove in L_i_toremove:
            dict_sample['L_contact_gimage'].pop(i_toremove)
            dict_sample['L_contact_ij_gimage'].pop(i_toremove)

        #Move grains (only Current)
        for grain in dict_sample['L_g']:
            if grain.group == 'Current' :
                grain.euler_semi_implicite(dict_algorithm['dt_DEM'])

        #periodic condition
        for grain in dict_sample['L_g']:
            #left wall
            if grain.center[0] < dict_sample['x_box_min'] :
                grain.center = grain.center.copy() + np.array([dict_sample['x_box_max'] - dict_sample['x_box_min'], 0])
                for i in range(len(grain.l_border)):
                    grain.l_border[i] = grain.l_border[i].copy() + np.array([dict_sample['x_box_max'] - dict_sample['x_box_min'], 0])
                    grain.l_border_x[i] = grain.l_border_x[i].copy() + dict_sample['x_box_max'] - dict_sample['x_box_min']
                #contact gimage needed to be convert into gg
                Owntools.convert_gimage_into_gg(grain, dict_sample, dict_material)
                #contact gg needed to be convert into gimage
                Owntools.convert_gg_into_gimage(grain, dict_sample, dict_material)
            #right wall
            elif grain.center[0] > dict_sample['x_box_max'] :
                grain.center = grain.center.copy() + np.array([dict_sample['x_box_min'] - dict_sample['x_box_max'], 0])
                for i in range(len(grain.l_border)):
                    grain.l_border[i] = grain.l_border[i].copy() + np.array([dict_sample['x_box_min'] - dict_sample['x_box_max'], 0])
                    grain.l_border_x[i] = grain.l_border_x[i].copy() + dict_sample['x_box_min'] - dict_sample['x_box_max']
                #contact gimage needed to be convert into gg
                Owntools.convert_gimage_into_gg(grain, dict_sample, dict_material)
                #contact gg needed to be convert into gimage
                Owntools.convert_gg_into_gimage(grain, dict_sample, dict_material)

        #Control the top group to have the pressure target
        dy_top, Fv = Control_Top_PID(dict_algorithm, dict_sollicitations['Vertical_Confinement_Force'], dict_sample['L_g'])
        #Apply confinement force
        for grain in dict_sample['L_g'] :
            if grain.group == 'Top':
                grain.move_as_a_group(np.array([0, dy_top]), dict_algorithm['dt_DEM'])
        dict_sample['y_box_max'] = dict_sample['y_box_max'] + dy_top

        #compute compacity, force applied on current grains and kinetic energy of current grains
        Surface_g = 0
        Force_applied = 0
        Ecin = 0
        for grain in dict_sample['L_g']:
            Surface_g = Surface_g + grain.surface
            if grain.group == 'Current':
                Force_applied = Force_applied + np.linalg.norm([grain.fx, grain.fy])
                Ecin = Ecin + 0.5 * grain.mass * np.dot(grain.v, grain.v)

        #tracker
        dict_tracker['compacity_L'].append(Surface_g/((dict_sample['y_box_max']-dict_sample['y_box_min'])*(dict_sample['x_box_max']-dict_sample['x_box_min'])))
        dict_tracker['vertical_force_L'].append(Fv)
        dict_tracker['dy_top_L'].append(dy_top)
        dict_tracker['Ecin_L'].append(Ecin)
        dict_tracker['Force_L'].append(Force_applied)

        if dict_algorithm['i_DEM'] % dict_algorithm['i_print_plot'] == 0:
            print('i_DEM',dict_algorithm['i_DEM'],': Confinement',int(100*Fv/dict_sollicitations['Vertical_Confinement_Force']),'%')
            if dict_algorithm['Debug_DEM'] :
                Owntools.Plot.Plot_Config_Confinement(dict_sample,dict_algorithm['i_DEM'])

        #Check stop conditions for DEM
        if dict_sample['L_g'] == []: #no more grains
            DEM_loop_statut = False
        if dict_algorithm['i_DEM'] >= dict_sollicitations['i_DEM_stop'] : #too many iterations
            DEM_loop_statut = False

    #plot trackers
    Owntools.Plot.Plot_own(list(range(0,len(dict_tracker['vertical_force_L']))),dict_tracker['vertical_force_L'], 'Vertical force', 'Debug/Confinement/vertical_force.png')
    Owntools.Plot.Plot_own(list(range(0,len(dict_tracker['dy_top_L']))),dict_tracker['dy_top_L'], 'dy', 'Debug/Confinement/dy.png')
    Owntools.Plot.Plot_own(list(range(0,len(dict_tracker['Ecin_L']))),dict_tracker['Ecin_L'], 'Kinetic energy', 'Debug/Confinement/ecin.png')
    Owntools.Plot.Plot_own(list(range(0,len(dict_tracker['Force_L']))),dict_tracker['Force_L'], 'Force applied', 'Debug/Confinement/force.png')

#-------------------------------------------------------------------------------

def Control_Top_PID(dict_algorithm, Force_target, L_g):
    """
    Control the upper wall to apply force.

    A PID corrector is applied.
        Input :
            an algorithm dictionnary (a dict)
            a confinement value (a float)
            a list of grain (a list)
        Output :
            the displacement of the top group (a float)
            a force applied on the top group before control (a float)
    """
    #compute vertical force applied on top group
    F = 0
    for grain in L_g :
        if grain.group == 'Top':
            F = F + grain.fy
    #compare with the target value
    error = F - Force_target #to have dy_top < 0 is F < Force_target
    #corrector
    ki = 0
    kd = 0
    dy_top = error * dict_algorithm['kp']
    #compare with maximum value
    if abs(dy_top) > dict_algorithm['dy_top_max'] :
        dy_top = np.sign(dy_top)*dict_algorithm['dy_top_max']

    return dy_top, F
