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

def Plot_Config_Sheared(dict_sample,i):
    """
    Plot loaded configuration during shearing.

        Input :
            a sample dictionnary (a dict)
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
    plt.savefig('Debug/Shear/Config_'+str(i)+'.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_Config_Confinement(dict_sample,i):
    """
    Plot loaded configuration during confinement.

        Input :
            a sample dictionnary (a dict)
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
    plt.savefig('Debug/Confinement/Config_'+str(i)+'.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_Contact(dict_sample,i):
    """
    Plot contact network.

        Input :
            a sample dictionnary (a dict)
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
                plt.plot(grain.l_border_x,grain.l_border_y, L_color_group[i_group])
    for grain in dict_sample['L_g_image']:
        for i_group in range(len(L_group)):
            if grain.group == L_group[i_group]:
                plt.plot(grain.l_border_x,grain.l_border_y,'-.', color = L_color_group[i_group])
    for contact in dict_sample['L_contact']+dict_sample['L_contact_gimage']:
        plt.plot([contact.g1.center[0], contact.g2.center[0]], [contact.g1.center[1], contact.g2.center[1]],'k')
    plt.axis('equal')
    plt.savefig('Debug/Shear/Contact_'+str(i)+'.png')
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
    plt.plot(dict_tracker['shear_L'], dict_tracker['vertical_force_before_L'])
    plt.plot([dict_tracker['shear_L'][0], dict_tracker['shear_L'][-1]], [dict_sollicitations['Vertical_Confinement_Force'], dict_sollicitations['Vertical_Confinement_Force']])
    plt.xlabel('Shear strain (-)')
    plt.ylabel('Vertical confinement force (µN)')

    #focus
    perc_confinement_L = []
    for force in dict_tracker['vertical_force_before_L'] :
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

#-------------------------------------------------------------------------------

def Plot_confinement_algorithm(dict_tracker):
    """
    Plot three plots to analyze the confinement algorithm.

    Plot :
        Evolution of the number of iterations needed to converge
        Evolution of the top group displacement
        Evolution of the difference between after and before vertical force

        Input :
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    plt.figure(1, figsize = (16,9))

    plt.subplot(131)
    plt.plot(dict_tracker['n_iteration_control_y_max_L'],'x')
    plt.xlabel('Iteration (-)')
    plt.ylabel('Number of iterations for the algorithm (-)')

    plt.subplot(132)
    plt.plot(dict_tracker['dy_top_L'],'x')
    plt.xlabel('Iteration (-)')
    plt.ylabel('Top group displacement (-)')

    plt.subplot(133)
    delta_force_L = []
    for i in range(len(dict_tracker['vertical_force_before_L'])):
        delta_force_L.append(dict_tracker['vertical_force_after_L'][i] - dict_tracker['vertical_force_before_L'][i])
    plt.plot(delta_force_L,'x')
    plt.xlabel('Iteration (-)')
    plt.ylabel('After - Before Fv (-)')

    plt.savefig('Debug/confinement_algorithm.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_own(x, y, title, path):
    """
    to write
    """
    plt.figure(1, figsize = (16,9))
    plt.plot(x,y)
    plt.title(title)
    plt.savefig(path)
    plt.close(1)
