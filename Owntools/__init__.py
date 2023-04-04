# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairies
#-------------------------------------------------------------------------------

from scipy.ndimage import binary_dilation
import numpy as np
import math

#Own
import Grain
import Contact_gg
import Contact_gimage
import Owntools.Plot

#-------------------------------------------------------------------------------
#functions
#-------------------------------------------------------------------------------

def convert_gimage_into_gg(grain, dict_sample, dict_material):
    """
    Convert a contact grain-image in a contact grain-grain.

        Input :
            a grain (a grain_tempo)
            a sample dictionnary (a dict)
            a material dictionnary (a dict)
        Output :
            Nothing, but the sample dictionnary is updated
    """
    L_i_contact_to_delete = []
    for ij_gimage in dict_sample['L_contact_ij_gimage'] :
        i_contact_gimage = dict_sample['L_contact_ij_gimage'].index(ij_gimage)
        if ij_gimage[0] == grain.id or ij_gimage[1] == grain.id:
            if ij_gimage[0] > ij_gimage[1] :
                ij_gg = (ij_gimage[1], ij_gimage[0])
            else :
                ij_gg = ij_gimage
            i_grain = 0
            grain_i = dict_sample['L_g'][i_grain]
            while not grain_i.id == ij_gg[0] :
                i_grain = i_grain + 1
                grain_i = dict_sample['L_g'][i_grain]
            j_grain = i_grain + 1
            grain_j = dict_sample['L_g'][j_grain]
            while not grain_j.id == ij_gg[1] :
                j_grain = j_grain + 1
                grain_j = dict_sample['L_g'][j_grain]
            if ij_gg not in dict_sample['L_contact_ij'] :
                #creation of contact
                dict_sample['L_contact_ij'].append(ij_gg)
                dict_sample['L_contact'].append(Contact_gg.Contact(dict_sample['id_contact'], grain_i, grain_j, dict_material))
                dict_sample['id_contact'] = dict_sample['id_contact'] + 1
                #transmit data
                dict_sample['L_contact'][-1].convert_gimage_in_gg(dict_sample['L_contact_gimage'][i_contact_gimage])
                #update neighborhood
                grain_i.neighborhood.append(grain_j)
            L_i_contact_to_delete.append(i_contact_gimage)
    #delete previous contact gimage
    L_i_contact_to_delete.reverse()
    for i_contact_to_delete in L_i_contact_to_delete :
        dict_sample['L_contact_gimage'].pop(i_contact_to_delete)
        dict_sample['L_contact_ij_gimage'].pop(i_contact_to_delete)

#-------------------------------------------------------------------------------

def convert_gg_into_gimage(grain, dict_sample, dict_material):
    """
    Convert a contact grain-grain in a contact grain-gimage.

        Input :
            a grain (a grain_tempo)
            a sample dictionnary (a dict)
            a material dictionnary (a dict)
        Output :
            Nothing, but the sample dictionnary is updated
    """
    L_i_contact_to_delete = []
    for ij_gg in dict_sample['L_contact_ij'] :
        i_contact_ij = dict_sample['L_contact_ij'].index(ij_gg)
        if ij_gg[0] == grain.id or ij_gg[1] == grain.id:
            ij_gimage = ij_gg
            if ij_gimage[1] in dict_sample['L_i_image'] :
                #contact gimage 1
                i_grain = 0
                grain = dict_sample['L_g'][i_grain]
                while not grain.id == ij_gimage[0] :
                    i_grain = i_grain + 1
                    grain = dict_sample['L_g'][i_grain]
                i_image = 0
                image = dict_sample['L_g_image'][i_image]
                while not image.id == ij_gimage[1] :
                    i_image = i_image + 1
                    image = dict_sample['L_g_image'][i_image]
                #creation of contact
                dict_sample['L_contact_ij_gimage'].append(ij_gimage)
                dict_sample['L_contact_gimage'].append(Contact_gimage.Contact_Image(dict_sample['id_contact'], grain, image, dict_material))
                dict_sample['id_contact'] = dict_sample['id_contact'] + 1
                #transmit data
                dict_sample['L_contact_gimage'][-1].convert_gg_in_gimage(dict_sample['L_contact'][i_contact_ij])
                #update neighborhood
                grain.neighborhood_image.append(image)

            #contact gimage 2
            ij_gimage = (ij_gg[1], ij_gg[0])
            if ij_gimage[1] in dict_sample['L_i_image'] :
                #contact gimage 1
                i_grain = 0
                grain = dict_sample['L_g'][i_grain]
                while not grain.id == ij_gimage[0] :
                    i_grain = i_grain + 1
                    grain = dict_sample['L_g'][i_grain]
                i_image = 0
                image = dict_sample['L_g_image'][i_image]
                while not image.id == ij_gimage[1] :
                    i_image = i_image + 1
                    image = dict_sample['L_g_image'][i_image]
                #creation of contact
                dict_sample['L_contact_ij_gimage'].append(ij_gimage)
                dict_sample['L_contact_gimage'].append(Contact_gimage.Contact_Image(dict_sample['id_contact'], grain, image, dict_material))
                dict_sample['id_contact'] = dict_sample['id_contact'] + 1
                #transmit data
                dict_sample['L_contact_gimage'][-1].convert_gg_in_gimage(dict_sample['L_contact'][i_contact_ij])
                #update neighborhood
                grain.neighborhood_image.append(image)

            L_i_contact_to_delete.append(i_contact_ij)
    #delete previous contact gimage
    L_i_contact_to_delete.reverse()
    for i_contact_to_delete in L_i_contact_to_delete :
        dict_sample['L_contact'].pop(i_contact_to_delete)
        dict_sample['L_contact_ij'].pop(i_contact_to_delete)

#-------------------------------------------------------------------------------

def convert_ic_to_real(dict_ic, dict_sample):
    """
    Convert the tempo grain from dict_ic into real grain in dict_sample.

        Input :
            an initial condition dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            the sample dictionnary gets grains information

    """
    dict_sample['L_g'] = []
    for grain_tempo in dict_ic['L_g_tempo'] :
        dict_sample['L_g'].append(Grain.Grain(grain_tempo))

#-------------------------------------------------------------------------------

def create_mesh(dict_algorithm, dict_material, dict_sample):
    """
    Create a mesh for the simulation.

    Material dictionnary is updated.

        Input :
            an algorithm dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            a mesh is generated in the dict_sample (a nx x ny numpy array)
            a material dictionnary is updated (two floats)
    """
    #create mesh and update dictionnary
    dict_sample['x_L'] = np.linspace(dict_sample['x_box_min'], dict_sample['x_box_max'], dict_algorithm['n_x'])
    dict_sample['y_L'] = np.linspace(dict_sample['y_box_min'], dict_sample['y_box_max'], dict_algorithm['n_y'])

    #plot mesh
    if dict_algorithm['Debug'] :
        Owntools.Plot.Plot_mesh(dict_sample)

    #From those date, add variables into material dict
    w = 4*math.sqrt((dict_sample['x_L'][1]-dict_sample['x_L'][0])**2+(dict_sample['y_L'][1]-dict_sample['y_L'][0])**2)
    double_well_height = 10*dict_material['kc_pf']/w/w
    dict_material['w'] = w
    dict_material['double_well_height'] = double_well_height

#-------------------------------------------------------------------------------

def Cosine_Profile(R,r,w):
    '''
    Compute the phase field variable at some point.

    A cosine profile is assumed (see https://mooseframework.inl.gov/source/ics/SmoothCircleIC.html).

    Input :
        the radius R of the grain in the direction (a float)
        the distance r between the current point and the center (a float)
        the width w of the interface (a float)
    Output :
        the value of the phase field variable (a float)
    '''
    #inside the grain
    if r<R-w/2:
        return 1
    #outside the grain
    elif r>R+w/2:
        return 0
    #inside the interface
    else :
        return 0.5*(1 + math.cos(math.pi*(r-R+w/2)/w))

#-------------------------------------------------------------------------------

def Interpolate_solute_out_grains(dict_algorithm, dict_sample) :
    """
    For all node of the mesh, if there is solute in one grain it is moved to the nearest node out of the grain.

    A node in multiple grain (contact area) is not considered in a grain a solute can exist.

        Input :
            an algorithm dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing, but the solute map is updated (a ny x nx numpy array)
    """
    #Create maps of different available nodes
    node_available_M = np.array(np.ones((len(dict_sample['y_L']),len(dict_sample['x_L'])),dtype = bool))
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            n_grain = 0
            for etai in dict_sample['L_etai']:
                if etai.etai_M[l][c] > 0.5 :
                    n_grain = n_grain + 1
            if n_grain == 1:
                node_available_M[l][c] = False

    #dilatation
    dilated_M = binary_dilation(node_available_M, dict_algorithm['struct_element'])

    #iteration on solute map to see if some solute is in not available position
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            if dict_sample['solute_M'][l][c] > 0 and not dilated_M[l][c]:
                solute_moved = False
                size_window = 1
                while not solute_moved :
                    i_window = 0
                    while not solute_moved and i_window <= size_window:
                        n_node_available = 0

                        #Look to move horizontaly and vertically
                        if i_window == 0 :
                            top_available = False
                            down_available = False
                            left_available = False
                            right_available = False
                            #to the top
                            if l - size_window > 0:
                                top_available = dilated_M[l-size_window][c]
                                if dilated_M[l-size_window][c] :
                                    n_node_available = n_node_available + 1
                            #to the down
                            if l + size_window < len(dict_sample['y_L']):
                                down_available = dilated_M[l+size_window][c]
                                if dilated_M[l+size_window][c] :
                                    n_node_available = n_node_available + 1
                            #to the left
                            if c - size_window > 0:
                                left_available = dilated_M[l][c-size_window]
                                if dilated_M[l][c-size_window] :
                                    n_node_available = n_node_available + 1
                            #to the right
                            if c + size_window < len(dict_sample['x_L']):
                                right_available = dilated_M[l][c+size_window]
                                if dilated_M[l][c+size_window] :
                                    n_node_available = n_node_available + 1

                            #move solute if et least one node is available
                            if n_node_available != 0 :
                                #to the top
                                if top_available:
                                    dict_sample['solute_M'][l-size_window][c] = dict_sample['solute_M'][l-size_window][c] + dict_sample['solute_M'][l][c]/n_node_available
                                #to the down
                                if down_available:
                                    dict_sample['solute_M'][l+size_window][c] = dict_sample['solute_M'][l+size_window][c] + dict_sample['solute_M'][l][c]/n_node_available
                                #to the left
                                if left_available:
                                    dict_sample['solute_M'][l][c-size_window] = dict_sample['solute_M'][l][c-size_window] + dict_sample['solute_M'][l][c]/n_node_available
                                #to the right
                                if right_available:
                                    dict_sample['solute_M'][l][c+size_window] = dict_sample['solute_M'][l][c+size_window] + dict_sample['solute_M'][l][c]/n_node_available
                                dict_sample['solute_M'][l][c] = 0
                                solute_moved = True

                        #Look to move diagonally
                        else :
                            top_min_available = False
                            top_max_available = False
                            down_min_available = False
                            down_max_available = False
                            left_min_available = False
                            left_max_available = False
                            right_min_available = False
                            right_max_available = False
                            #to the top
                            if l - size_window > 0:
                                if c - i_window > 0 :
                                    top_min_available = dilated_M[l-size_window][c-i_window]
                                    if dilated_M[l-size_window][c-i_window] :
                                        n_node_available = n_node_available + 1
                                if c + i_window < len(dict_sample['x_L']):
                                    top_max_available = dilated_M[l-size_window][c+i_window]
                                    if dilated_M[l-size_window][c+i_window] :
                                        n_node_available = n_node_available + 1
                            #to the down
                            if l + size_window < len(dict_sample['y_L']):
                                if c - i_window > 0 :
                                    down_min_available = dilated_M[l+size_window][c-i_window]
                                    if dilated_M[l+size_window][c-i_window] :
                                        n_node_available = n_node_available + 1
                                if c + i_window < len(dict_sample['x_L']):
                                    down_max_available = dilated_M[l+size_window][c+i_window]
                                    if dilated_M[l+size_window][c+i_window] :
                                        n_node_available = n_node_available + 1
                            #to the left
                            if c - size_window > 0:
                                if l - i_window > 0 :
                                    left_min_available = dilated_M[l-i_window][c-size_window]
                                    if dilated_M[l-i_window][c-size_window] :
                                        n_node_available = n_node_available + 1
                                if l + i_window < len(dict_sample['y_L']):
                                    left_max_available = dilated_M[l+i_window][c-size_window]
                                    if dilated_M[l+i_window][c-size_window] :
                                        n_node_available = n_node_available + 1
                            #to the right
                            if c + size_window < len(dict_sample['x_L']):
                                if l - i_window > 0 :
                                    right_min_available = dilated_M[l-i_window][c+size_window]
                                    if dilated_M[l-i_window][c+size_window] :
                                        n_node_available = n_node_available + 1
                                if l + i_window < len(dict_sample['y_L']):
                                    right_max_available = dilated_M[l+i_window][c+size_window]
                                    if dilated_M[l+i_window][c+size_window] :
                                        n_node_available = n_node_available + 1

                            #move solute if at least one node is available
                            if n_node_available != 0 :
                                #to the top
                                if top_min_available:
                                    dict_sample['solute_M'][l-size_window][c-i_window] = dict_sample['solute_M'][l-size_window][c-i_window] + dict_sample['solute_M'][l][c]/n_node_available
                                if top_max_available:
                                    dict_sample['solute_M'][l-size_window][c+i_window] = dict_sample['solute_M'][l-size_window][c+i_window] + dict_sample['solute_M'][l][c]/n_node_available
                                #to the down
                                if down_min_available:
                                    dict_sample['solute_M'][l+size_window][c-i_window] = dict_sample['solute_M'][l+size_window][c-i_window] + dict_sample['solute_M'][l][c]/n_node_available
                                if down_max_available:
                                    dict_sample['solute_M'][l+size_window][c+i_window] = dict_sample['solute_M'][l+size_window][c+i_window] + dict_sample['solute_M'][l][c]/n_node_available
                                #to the left
                                if left_min_available:
                                    dict_sample['solute_M'][l-i_window][c-size_window] = dict_sample['solute_M'][l-i_window][c-size_window] + dict_sample['solute_M'][l][c]/n_node_available
                                if left_max_available:
                                    dict_sample['solute_M'][l+i_window][c-size_window] = dict_sample['solute_M'][l+i_window][c-size_window] + dict_sample['solute_M'][l][c]/n_node_available
                                #to the right
                                if right_min_available:
                                    dict_sample['solute_M'][l-i_window][c+size_window] = dict_sample['solute_M'][l-i_window][c+size_window] + dict_sample['solute_M'][l][c]/n_node_available
                                if right_max_available:
                                    dict_sample['solute_M'][l+i_window][c+size_window] = dict_sample['solute_M'][l+i_window][c+size_window] + dict_sample['solute_M'][l][c]/n_node_available
                                dict_sample['solute_M'][l][c] = 0
                                solute_moved = True

                        i_window = i_window + 1
                    size_window = size_window + 1
