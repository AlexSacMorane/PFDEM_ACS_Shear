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
import Owntools.Plot

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def DEM_shear_load(dict_algorithm, dict_material, dict_sample, dict_sollicitation, simulation_report):
    """
    Loading the granular system with vertical load and shear.

        Input :
            an algorithm dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitation dictionnary (a dict)
            a simultion report (a report)
        Output :
            Nothing, but sample dictionnary is updated
    """
    i_DEM_0 = dict_sample['i_DEM']
    Shear_strain = 0
    #compute the sample height
    min_value = min(dict_sample['L_g'][0].l_border_y)
    max_value = max(dict_sample['L_g'][0].l_border_y)
    for grain in dict_sample['L_g']:
        if min(grain.l_border_y) < min_value :
            min_value = min(grain.l_border_y)
        if max(grain.l_border_y) > max_value :
            max_value = max(grain.l_border_y)
    Sample_height = max_value - min_value
    #track total displacement of grains
    for grain in dict_sample['L_g'] :
        grain.track_u = True
        grain.total_ux = 0
        grain.total_uy = 0
    DEM_loop_statut = True
    #Initialisation
    dict_sample['L_contact'] = []
    dict_sample['L_contact_ij'] = []
    dict_sample['L_contact_gimage'] = []
    dict_sample['L_contact_ij_gimage'] = []
    dict_sample['id_contact'] = 0
    dict_sample['L_g_image'] = []
    dict_sample['L_i_image'] = []

    while DEM_loop_statut :

        dict_sample['i_DEM'] = dict_sample['i_DEM'] + 1

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
        if (dict_sample['i_DEM']-i_DEM_0-1) % dict_sample['i_update_neighborhoods_com']  == 0:
            Contact_gg.Update_Neighborhoods(dict_sample)
            Contact_gimage.Update_Neighborhoods(dict_sample)
        Contact_gg.Grains_contact_Neighborhoods(dict_sample,dict_material)
        Contact_gimage.Grains_contact_Neighborhoods(dict_sample,dict_material)

        #Sollicitation computation
        for grain in dict_sample['L_g']:
             grain.init_F_control(dict_sollicitation['gravity'])
        for contact in  dict_sample['L_contact']+dict_sample['L_contact_gimage']+dict_sample['L_contact_gw']:
            contact.normal()
            contact.tangential(dict_sample['dt_DEM_IC'])

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
                grain.euler_semi_implicite(dict_sample['dt_DEM_IC'],10*dict_sample['Ecin_ratio_IC'])

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
        dy_top, Fv = Control_Top_NR(dict_sollicitation['Vertical_Confinement_Force'],dict_sample['L_contact']+dict_sample['L_contact_gimage'],dict_sample['L_g'])
        dict_sample['y_box_max'] = dict_sample['y_box_max'] + dy_top
        for grain in dict_sample['L_g'] :
            if grain.group == 'Top':
                grain.move_as_a_group(np.array([dict_sollicitation['Shear_velocity']*dict_sample['dt_DEM_IC'], dy_top]), dict_sample['dt_DEM_IC'])
        #Update shear strain
        Shear_strain = Shear_strain + dict_sollicitation['Shear_velocity']*dict_sample['dt_DEM_IC'] / Sample_height

        if dict_sample['i_DEM'] % dict_sample['i_print_plot'] == 0:
            print('i_DEM',dict_sample['i_DEM'],'and Confinement',int(100*Fv/dict_sollicitation['Vertical_Confinement_Force']),'% and Shear',int(100*Shear_strain/dict_sollicitation['Shear_strain_target']),'%')
            if dict_sample['Debug_DEM'] :
                Owntools.Plot_Config_Loaded(dict_sample,dict_sample['i_DEM'])

        #Check stop conditions for DEM
        if Shear_strain >= dict_sollicitation['Shear_strain_target'] :
             DEM_loop_statut = False
        if dict_sample['L_g'] == []:
            DEM_loop_statut = False

    #plot total displacement field
    Owntools.Plot_total_U(dict_sample)

#-------------------------------------------------------------------------------

def Control_Top_NR(Force_target,L_contact_gg,L_g):
    """
    Control the upper wall to apply force.

    A Newton-Raphson method is applied to verify the confinement.
        Input :
            a confinement value (a float)
            a list of contact grain - grain and grain - image (a list)
            a list of temporary grain (a list)
        Output :
            the displacement of the top group (a float)
            a force applied on the top group before control (a float)
    """
    F = 0
    overlap_L = []
    k_L = []
    for contact in L_contact_gg:
        #compute force applied, save contact overlap and spring
        #only the normal component is considered, no tangential one, no damping
        if (contact.g1.group == 'Current' and contact.g2.group == 'Top') :
            F = F - contact.F_2_1_n * contact.pc_normal[1] #- because F_2_1_n < 0
            overlap_L.append(contact.overlap_normal)
            k_L.append(contact.k*contact.pc_normal[1])
        elif (contact.g1.group == 'Top' and contact.g2.group == 'Current') :
            F = F + contact.F_2_1_n * contact.pc_normal[1] #+ because F_2_1_n < 0 and pc_normal from top to current
            overlap_L.append(contact.overlap_normal)
            k_L.append(-contact.k*contact.pc_normal[1]) #because pc_normal from top to current

    if overlap_L != []:
        i_NR = 0
        dy_top = 0
        ite_criteria = True
        #control the upper wall
        if -0.01*Force_target<error_on_ymax_f(dy_top,overlap_L,k_L,Force_target) and error_on_ymax_f(dy_top,overlap_L,k_L,Force_target)<0.01*Force_target:
            ite_criteria = False
        while ite_criteria :
            i_NR = i_NR + 1
            dy_top = dy_top - error_on_ymax_f(dy_top,overlap_L,k_L,Force_target)/error_on_ymax_df(dy_top,overlap_L,k_L)
            if i_NR > 100: #Maximum try
                ite_criteria = False
            if -0.01*Force_target<error_on_ymax_f(dy_top,overlap_L,k_L,Force_target) and error_on_ymax_f(dy_top,overlap_L,k_L,Force_target)<0.01*Force_target:
                ite_criteria = False

    else :
        #if there is no contact with the upper wall, the wall is reset
        dy_top = Reset_y_max(L_g)

    return dy_top, F

#-------------------------------------------------------------------------------

def error_on_ymax_f(dy,overlap_L,k_L,Force_target) :
    """
    Compute the function f to control the upper wall. It is the difference between the force applied and the target value.

        Input :
            an increment of the upper wall position (a float)
            a list of overlap for contact between temporary grain and upper wall (a list)
            a list of spring for contact between temporary grain and upper wall (a list)
            a confinement force (a float)
        Output :
            the difference between the force applied and the confinement (a float)
    """
    f = Force_target
    for i in range(len(overlap_L)):
        f = f - k_L[i]*(max(overlap_L[i]-dy,0))**(3/2)
    return f

#-------------------------------------------------------------------------------

def error_on_ymax_df(dy,overlap_L,k_L) :
    """
    Compute the derivative function df to control the upper wall (error_on_ymax_f()).

        Input :
            an increment of the upper wall position (a float)
            a list of overlap for contact between temporary grain and upper wall (a list)
            a list of spring for contact between temporary grain and upper wall (a list)
        Output :
            the derivative of error_on_ymax_f() (a float)
    """
    df = 0
    for i in range(len(overlap_L)):
        df = df + 3/2*k_L[i]*(max(overlap_L[i]-dy,0))**(1/2)
    return df

#-------------------------------------------------------------------------------

def Reset_y_max(L_g):
    """
    The Top group is going down by a fraction of the mean radius.

    The confinement force is not verified.

        Input :
            the list of temporary grains (a list)
        Output :
            the upper wall position (a float)
    """
    R_mean = 0
    n_grain = 0
    for grain in L_g :
        if grain.group == 'Top' :
            R_mean = R_mean + grain.radius
            n_grain = n_grain + 1
    dy_top = - R_mean / n_grain / 20
    return dy_top
