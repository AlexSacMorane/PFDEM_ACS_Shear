# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions used to compute parameters or variables in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np
import math
import random
from scipy.ndimage import binary_dilation

#-------------------------------------------------------------------------------

def Compute_Emec(dict_material, dict_sample, dict_sollicitations):
    '''
    Compute the mechanical energy over the sample.

    There is an iteration over all the contacts detected (grain-grain). First, the sum of the minimal grain phase variables is computed.
    Then, the mechanical energy is computed.

        Input :
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitation dictionnary (a dict)
        Output :
            Nothing, but the sample dictionnary gets an updated value for the mechanical term (a nx x ny numpy array)
    '''
    #Initialisation
    Emec_M = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))))
    #contact grain-grain part
    for contact in dict_sample['L_contact']:
        #extract a spatial zone
        x_min = min(min(contact.g1.l_border_x),min(contact.g2.l_border_x))-dict_material['w']
        x_max = max(max(contact.g1.l_border_x),max(contact.g2.l_border_x))+dict_material['w']
        y_min = min(min(contact.g1.l_border_y),min(contact.g2.l_border_y))-dict_material['w']
        y_max = max(max(contact.g1.l_border_y),max(contact.g2.l_border_y))+dict_material['w']
        #look for this part inside the global mesh
        #create search list
        x_L_search_min = abs(np.array(dict_sample['x_L'])-x_min)
        x_L_search_max = abs(np.array(dict_sample['x_L'])-x_max)
        y_L_search_min = abs(np.array(dict_sample['y_L'])-y_min)
        y_L_search_max = abs(np.array(dict_sample['y_L'])-y_max)
        #get index
        i_x_min = list(x_L_search_min).index(min(x_L_search_min))
        i_x_max = list(x_L_search_max).index(min(x_L_search_max))
        i_y_min = list(y_L_search_min).index(min(y_L_search_min))
        i_y_max = list(y_L_search_max).index(min(y_L_search_max))
        #Initialisation
        sum_min_etai = 0
        #compute the sum over the sample of the minimum of etai
        for l in range(i_y_min, i_y_max):
            for c in range(i_x_min, i_x_max):
                sum_min_etai = sum_min_etai + min(contact.g1.etai_M[-1-l][c],contact.g2.etai_M[-1-l][c])
        if sum_min_etai != 0 :
            #compute the variable e_mec
            e_mec = dict_sollicitations['alpha']/sum_min_etai
        else :
            e_mec = 0
        #compute the distribution of the mechanical energy
        for l in range(i_y_min, i_y_max):
            for c in range(i_x_min, i_x_max):
                Emec_M[-1-l][c] = Emec_M[-1-l][c] + e_mec*min(contact.g1.etai_M[-1-l][c],contact.g2.etai_M[-1-l][c])

    #Update element in dictionnary
    dict_sample['Emec_M'] = Emec_M

#-------------------------------------------------------------------------------

def Compute_kc(dict_algorithm, dict_material, dict_sample):
    '''
    Compute the solute diffusion coefficient field in the sample.

    Here, a dilation method is applied. For all node, a Boolean variable is defined.
    This variable is True at least 2 eta_i are greater than 0.5 (in the contact zone).
                  is True all etai are lower than 0.5 (in the pore zone).
                  is False else.

    A dilation method is applied, the size of the structural element is the main case.

    The diffusion map is built on the Boolean map. If the variable is True, the diffusion is kc, else 0.

        Input :
            an algorithm dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing but the dictionnary gets an updated value for the solute diffusion coefficient map (a nx x ny numpy array)
    '''
    #Initialisation
    on_off_M = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))), dtype = bool)

    #compute the on off map
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            n_etai_over_0_5 = 0 #counter
            for etai in dict_sample['L_etai'] :
                if etai.etai_M[-1-l][c] > 0.5 :
                    n_etai_over_0_5 = n_etai_over_0_5 + 1
            #at a contact
            if n_etai_over_0_5 >= 2:
                on_off_M[-l-1][c] = True
            #in the pore space
            elif n_etai_over_0_5 == 0:
                on_off_M[-l-1][c] = True

    #dilatation
    dilated_M = binary_dilation(on_off_M, dict_algorithm['struct_element'])

    #compute the map of the solute diffusion coefficient
    kc_M = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))))
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            if dilated_M[-1-l][c] :
                kc_M[-1-l][c] = dict_material['kappa_c']

    #Update element in dictionnary
    dict_sample['kc_M'] = kc_M
