# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairies
#-------------------------------------------------------------------------------

from pathlib import Path
import matplotlib.pyplot as plt
import imageio

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

def Plot_mesh(dict_sample):
    """
    Plot group distribution and mesh.

        Input :
            a sample disctionnary (a dict)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    L_color_group = ['k','r','b']
    L_group = ['Current', 'Bottom', 'Top']
    fig = plt.figure(1,figsize=(16,9.12))
    #x-axis
    for x in dict_sample['x_L']:
        plt.plot([x, x],[dict_sample['y_L'][0], dict_sample['y_L'][-1]],'k')
    #y-axis
    for y in dict_sample['y_L']:
        plt.plot([dict_sample['x_L'][0], dict_sample['x_L'][-1]],[y, y],'k')
    #grains
    for grain in dict_sample['L_g']:
        for i_group in range(len(L_group)):
            if grain.group == L_group[i_group] :
                plt.plot(grain.l_border_x, grain.l_border_y, L_color_group[i_group])
    plt.axis("equal")
    fig.savefig('Debug/Mesh.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_etais(dict_sample):
    """
    Plot etais distribution.

        Input :
            a sample disctionnary (a dict)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    fig = plt.figure(1,figsize=(16,9.12))

    L_color = ['red', 'royalblue', 'forestgreen', 'gold', 'hotpink', 'skyblue', 'chocolate', 'darkkhaki', 'darkorchid', 'silver']
    #etas
    etas_M = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))))
    for etai in dict_sample['L_etai']:
        etas_M = etas_M + etai.etai_M
    im = plt.imshow(etas_M, interpolation='nearest', extent=[min(dict_sample['x_L']),max(dict_sample['x_L']), min(dict_sample['y_L']),max(dict_sample['y_L'])], vmin = 0, vmax = 1)
    plt.colorbar(im)
    #grains
    for grain in dict_sample['L_g']:
        plt.plot(grain.l_border_x,grain.l_border_y,color=L_color[grain.etai],linewidth=3)
    plt.title('Phase field and grains')
    plt.axis('equal')
    plt.xlim(min(dict_sample['x_L']),max(dict_sample['x_L']))
    plt.axis("equal")
    fig.savefig('Debug/Etais.png')
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
    #solute
    im = plt.imshow(dict_sample['solute_M'], interpolation='nearest', extent=[min(dict_sample['x_L']),max(dict_sample['x_L']), min(dict_sample['y_L']),max(dict_sample['y_L'])], vmin = 0, vmax = 0.1)
    plt.colorbar(im)
    #grains
    for grain in dict_sample['L_g']:
        for i_group in range(len(L_group)):
            if grain.group == L_group[i_group] :
                plt.plot(grain.l_border_x,grain.l_border_y,L_color_group[i_group])
    for grain in dict_sample['L_g_image']:
        for i_group in range(len(L_group)):
            if grain.group == L_group[i_group]:
                plt.plot(grain.l_border_x,grain.l_border_y,'-.', color = L_color_group[i_group])
    plt.title('Solute and grains')
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

def Plot_mp4(template_name,name_video):
    '''The goal of this function is to create a movie with pictures.

    from https://www.blog.pythonlibrary.org/2021/06/23/creating-an-animated-gif-with-python/

        Input :
            the template of the pictures used (a string)
            the name of the video (a string)
        Output :
            a movie file (a .mp4)
    '''
    #look for the largest iteration
    i_f = 0
    plotpath = Path(template_name+str(i_f)+'.png')
    while plotpath.exists():
        i_f = i_f + 1
        plotpath = Path(template_name+str(i_f)+'.png')

    fileList = []
    for i in range(0,i_f):
        fileList.append(template_name+str(i)+'.png')

    duration_movie  = 10 #sec
    writer = imageio.get_writer(name_video, fps=int(i_f/duration_movie))
    for im in fileList:
        writer.append_data(imageio.imread(im))
    writer.close()

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

def Plot_strain_dy_top(dict_algorithm, dict_tracker):
    """
    Plot the curve strain - dy of the top group.

        Input :
            an algorithm dictionnary (a dict)
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    plt.figure(1, figsize = (16,9))
    plt.plot(dict_tracker['shear_L'], dict_tracker['dy_top_L'])
    plt.plot([dict_tracker['shear_L'][0],dict_tracker['shear_L'][-1]], [-dict_algorithm['dy_top_max'], -dict_algorithm['dy_top_max']], 'r')
    plt.plot([dict_tracker['shear_L'][0],dict_tracker['shear_L'][-1]], [ dict_algorithm['dy_top_max'],  dict_algorithm['dy_top_max']], 'r')
    plt.xlabel('Shear strain (-)')
    plt.ylabel('dy Top group (µm)')
    plt.savefig('Debug/strain_dy_top.png')
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
