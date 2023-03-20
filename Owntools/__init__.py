# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairies
#-------------------------------------------------------------------------------

import Grain
import Contact_gg
import Contact_gimage

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
