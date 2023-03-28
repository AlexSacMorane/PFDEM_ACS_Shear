# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the algorithm to apply confinement force.
"""

#-------------------------------------------------------------------------------
#Librairies
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

#Own
import Grain
import Contact_gg
import Contact_gimage

#-------------------------------------------------------------------------------
#Classes
#-------------------------------------------------------------------------------

class Grain_copy(Grain.Grain):
    """
    A copy of a polygnal real grain used in the simulation.
    """

#-------------------------------------------------------------------------------

class Grain_Image_copy(Grain.Grain_Image):
    """
    A copy of a polygnal image grain used in the simulation.
    """

#-------------------------------------------------------------------------------

class Contact_copy(Contact_gg.Contact):
    """
    A copy of a contact grain-grain used in the simulation.
    """

#-------------------------------------------------------------------------------

class Contact_Image_copy(Contact_gimage.Contact_Image):
    """
    A copy of a contact grain-image used in the simulation.
    """

#-------------------------------------------------------------------------------
#functions
#-------------------------------------------------------------------------------

def Control_y_max_copy(dict_algorithm, dict_material, dict_sample, dict_sollicitations):
    """
    Apply the confinement force.

    The top group and upper grain from current group are extracted.
    The top group is translated until the confinement force (sum of the component y of the force applied on grains) reachs the target value

        Input :
            an algorithm dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
        Output :
            sum_dy_top, fy_sum, n_iteration
    """
    #sort grains
    L_top_copy = []
    L_current_copy = []
    for grain in dict_sample['L_g']:
        if grain.group =='Top':
            L_top_copy.append(Grain_copy(grain))
        elif grain.group == 'Current' :
            if dict_sample['y_box_max'] - dict_algorithm['height_to_consider'] < grain.center[1] :
                L_current_copy.append(Grain_copy(grain))
    #find the new position of the top group to apply confinement
    conv = False
    sum_dy_top = 0
    n_switch_sign = 0
    Delta_fy = 0
    n_iteration = 0
    while not conv :
        n_iteration = n_iteration + 1
        #create image
        L_g_image = []
        for grain in L_top_copy + L_current_copy:
            #left wall
            if (grain.center[0] - dict_sample['x_box_min']) < dict_algorithm['d_to_image'] :
                L_g_image.append(Grain_Image_copy(grain, 'left'))
            #right wall
            elif (dict_sample['x_box_max'] - grain.center[0]) < dict_algorithm['d_to_image'] :
                L_g_image.append(Grain_Image_copy(grain, 'right'))
        #translate image
        for image in L_g_image:
            if image.position == 'left' :
                image.translation(np.array([dict_sample['x_box_max'] - dict_sample['x_box_min'], 0]))
            elif image.position == 'right' :
                image.translation(np.array([dict_sample['x_box_min'] - dict_sample['x_box_max'], 0]))

        i_contact = 0
        #contact grain-grain
        L_contact_gg = []
        for grain_current in L_current_copy :
            for grain_top in L_top_copy :
                if Contact_gg.Grains_Polyhedral_contact_f(grain_current, grain_top) :
                    L_contact_gg.append(Contact_copy(i_contact, grain_current, grain_top, dict_material))
                    i_contact = i_contact + 1
        for grain_current in L_current_copy :
            for grain_top in L_top_copy :
                if Contact_gg.Grains_Polyhedral_contact_f(grain_current, grain_top) :
                    L_contact_gg.append(Contact_copy(i_contact, grain_current, grain_top, dict_material))
                    i_contact = i_contact + 1

        #contact grain-image
        L_contact_gimage = []
        for grain_top in L_top_copy :
            for image in L_g_image :
                if image.group == 'Current' :
                    if Contact_gg.Grains_Polyhedral_contact_f(grain_top, image) :
                        L_contact_gimage.append(Contact_Image_copy(i_contact, grain_top, image, dict_material))
                        #the contacts grain(current)-image(top) are not considered as the contact apply a force only on the grain.
                        #A focus is done on the top group.
                        i_contact = i_contact + 1

        #Sollicitation computation
        for grain in L_top_copy:
             grain.init_F_control(0)
        for contact in L_contact_gg + L_contact_gimage:
            #do not consider the contact inside top and bottom groups
            if not (contact.g1.group == 'Top' and contact.g2.group =='Top') :
                contact.normal()
                contact.tangential(dict_algorithm['dt_DEM'])

        #compute the force applied on top group
        fy_sum = 0
        for grain in L_top_copy:
            fy_sum = fy_sum + grain.fy
        print(fy_sum)
        #compare with target value
        Delta_fy_before = Delta_fy
        Delta_fy = fy_sum - dict_sollicitations['Vertical_Confinement_Force']
        if Delta_fy*Delta_fy_before < 0 :
            n_switch_sign = n_switch_sign + 1
        if abs(Delta_fy) < 0.05*dict_sollicitations['Vertical_Confinement_Force']:
            conv = True
        else :
            dy_top = Delta_fy/abs(Delta_fy)*dict_algorithm['dy_top']/(2**n_switch_sign)
            sum_dy_top = sum_dy_top + dy_top
            #move top grains
            for grain in L_top_copy :
                grain.move_as_a_group(np.array([0, dy_top]), 1) #dt_dEM is set to 1, not used for top group

    return sum_dy_top, fy_sum, n_iteration
