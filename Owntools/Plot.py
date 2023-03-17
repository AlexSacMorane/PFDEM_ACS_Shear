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
    fig.savefig('Debug/Configuration/Group_Distribution.png')
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
    plt.savefig('Debug/Configuration/Shear/Config_Loaded_'+str(i)+'.png')
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
    plt.savefig('Debug/Configuration/Total_U_Sheared.png')
    plt.close(1)

#-------------------------------------------------------------------------------
