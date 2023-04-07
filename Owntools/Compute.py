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
