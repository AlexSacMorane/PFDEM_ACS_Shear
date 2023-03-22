# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairies
#-------------------------------------------------------------------------------

import matplotlib.pyplot as plt

#Own
import Grain
import Contact_gg
import Contact_gimage

#-------------------------------------------------------------------------------
#functions
#-------------------------------------------------------------------------------

def Plot_group_distribution(dict_ic):
    """
    Plot group distribution.

        Input :
            an initial configuration disctionnary (a dict)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    L_color_group = ['k','r','b']
    L_group = ['Current', 'Bottom', 'Top']
    fig = plt.figure(1,figsize=(16,9.12))
    for grain in dict_ic['L_g_tempo']:
        for i_group in range(len(L_group)):
            if grain.group == L_group[i_group] :
                plt.plot(grain.l_border_x, grain.l_border_y, L_color_group[i_group])
    plt.axis("equal")
    fig.savefig('Debug/Group_Distribution.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_Config_Loaded(dict_sample,i):
    """
    Plot loaded configuration.

        Input :
            a list of temporary grain (a list)
            the coordinates of the walls (four floats)
            an iteration (a int)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    L_color_group = ['k','r','b']
    L_group = ['Current', 'Bottom', 'Top']
    plt.figure(1,figsize=(16,9))
    for grain in dict_sample['L_g']:
        for i_group in range(len(L_group)):
            if grain.group == L_group[i_group] :
                plt.plot(grain.l_border_x,grain.l_border_y,L_color_group[i_group])
    for grain in dict_sample['L_g_image']:
        for i_group in range(len(L_group)):
            if grain.group == L_group[i_group]:
                plt.plot(grain.l_border_x,grain.l_border_y,'-.', color = L_color_group[i_group])
    plt.axis('equal')
    plt.savefig('Debug/Shear/Config_Loaded_'+str(i)+'.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_total_U(dict_sample):
    """
    Plot the map of the total displacement tracked.

        Input :
            a list of temporary grain (a list)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    L_color_group = ['k','r','b']
    L_group = ['Current', 'Bottom', 'Top']
    plt.figure(1,figsize=(16,9))
    x_L = []
    y_L = []
    u_L = []
    v_L = []
    for grain in dict_sample['L_g']:
        for i_group in range(len(L_group)):
            if grain.group == L_group[i_group] :
                plt.plot(grain.l_border_x,grain.l_border_y,L_color_group[i_group])
        if grain.group == 'Current' :
            x_L.append(grain.center[0])
            y_L.append(grain.center[1])
            u_L.append(grain.total_ux)
            v_L.append(grain.total_uy)
    plt.quiver(x_L,y_L,u_L,v_L)
    plt.axis('equal')
    plt.savefig('Debug/Total_U_Sheared.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_strain_compacity(dict_tracker):
    """
    Plot the curve strain - compacity.

        Input :
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    plt.figure(1, figsize = (16,9))
    plt.plot(dict_tracker['shear_L'], dict_tracker['compacity_L'])
    plt.xlabel('Shear strain (-)')
    plt.ylabel('Compacity (-)')
    plt.savefig('Debug/strain_compacity.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_strain_confinement(dict_tracker, dict_sollicitations):
    """
    Plot the curve strain - vertical force (confinement).

        Input :
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    plt.figure(1, figsize = (16,9))

    plt.subplot(121)
    plt.plot(dict_tracker['shear_L'], dict_tracker['vertical_force_L'])
    plt.plot([dict_tracker['shear_L'][0], dict_tracker['shear_L'][-1]], [dict_sollicitations['Vertical_Confinement_Force'], dict_sollicitations['Vertical_Confinement_Force']])
    plt.xlabel('Shear strain (-)')
    plt.ylabel('Vertical confinement force (µN)')

    #focus
    perc_confinement_L = []
    for force in dict_tracker['vertical_force_L'] :
        perc_confinement_L.append(force/dict_sollicitations['Vertical_Confinement_Force']*100)
    plt.subplot(122)
    plt.plot(dict_tracker['shear_L'], perc_confinement_L)
    plt.ylim((95,105))
    plt.xlabel('Shear strain (-)')
    plt.ylabel('Vertical confinement force (% of the target value)')

    plt.savefig('Debug/strain_confinement.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_strain_mu_sample(dict_tracker):
    """
    Plot the curve strain - sample friction.

        Input :
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    plt.figure(1, figsize = (16,9))
    plt.plot(dict_tracker['shear_L'], dict_tracker['mu_sample_L'])
    plt.xlabel('Shear strain (-)')
    plt.ylabel('Sample friction (-)')
    plt.savefig('Debug/strain_sample_friction.png')
    plt.close(1)
