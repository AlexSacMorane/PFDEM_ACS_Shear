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

def DEM_shear_load(dict_algorithm, dict_geometry, dict_material, dict_sample, dict_sollicitations, dict_tracker, simulation_report):
    """
    Loading the granular system with vertical load and shear.

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitation dictionnary (a dict)
            a tracker dictionnary (a dict)
            a simultion report (a report)
        Output :
            Nothing, but sample dictionnary is updated
    """
    dict_algorithm['i_DEM'] = 0
    i_DEM_0 = dict_algorithm['i_DEM']
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
    #compute the inertial number
    dict_sample['I_number'] = dict_sollicitations['Shear_velocity']/Sample_height*2*dict_geometry['R_mean']*math.sqrt(dict_material['rho_surf']*dict_sollicitations['Vertical_Confinement_Linear_Force'])
    simulation_report.write_and_print('Inertial number : '+str(dict_sample['I_number'])+'\n','Inertial number : '+str(dict_sample['I_number']))
    #must be under 10-3 to consider critical state
    simulation_report.write_and_print('Expected number of iterations : '+str(int(dict_sollicitations['Shear_strain_target']*Sample_height/(dict_sollicitations['Shear_velocity']*dict_algorithm['dt_DEM'])))+'\n\n','Expected number of iterations : '+str(int(dict_sollicitations['Shear_strain_target']*Sample_height/(dict_sollicitations['Shear_velocity']*dict_algorithm['dt_DEM'])))+'\n')
    #Switch on tracks
    for grain in dict_sample['L_g'] :
        #track rigid body motion to move pf
        grain.track_rbm_pf_interpolation = True
        grain.u_pf_interpolation = np.array([0,0])
        grain.dtheta_pf_interpolation = 0
        #track total displacement of grains (plot)
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
    #trackers
    dict_tracker['vertical_force_before_L'] = []
    dict_tracker['vertical_force_after_L'] = []
    #dict_tracker['n_iteration_control_y_max_L'] = []
    dict_tracker['dy_top_L'] = []
    dict_tracker['shear_L'] = []
    dict_tracker['mu_L'] = []
    dict_tracker['compacity_L'] = []
    dict_tracker['mu_sample_L'] = []

    while DEM_loop_statut :

        dict_algorithm['i_DEM'] = dict_algorithm['i_DEM'] + 1

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
        if (dict_algorithm['i_DEM']-i_DEM_0-1) % dict_algorithm['i_update_neighborhoods']  == 0:
            Contact_gg.Update_Neighborhoods(dict_algorithm, dict_sample)
            Contact_gimage.Update_Neighborhoods(dict_algorithm, dict_sample)
        Contact_gg.Grains_contact_Neighborhoods(dict_sample,dict_material)
        Contact_gimage.Grains_contact_Neighborhoods(dict_sample,dict_material)

        #Sollicitation computation
        for grain in dict_sample['L_g']:
             grain.init_F_control(dict_sollicitations['gravity'])
        for contact in dict_sample['L_contact']+dict_sample['L_contact_gimage']:
            #do not consider the contact inside top and bottom groups
            if not (contact.g1.group == 'Top' and contact.g2.group =='Top') and not (contact.g1.group == 'Bottom' and contact.g2.group =='Bottom') :
                contact.normal()
                contact.tangential(dict_algorithm['dt_DEM'])

        #Move grains (only Current)
        for grain in dict_sample['L_g']:
            if grain.group == 'Current' :
                grain.euler_semi_implicite(dict_algorithm['dt_DEM'])

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

        #Compute the sample friction coefficient
        sum_fx_top = 0
        sum_fy_top = 0
        for grain in dict_sample['L_g']:
            if grain.group == 'Top':
                sum_fx_top = sum_fx_top + grain.fx
                sum_fy_top = sum_fy_top + grain.fy
        if sum_fy_top != 0 : #else keep same value
            mu_sample = abs(sum_fx_top / sum_fy_top)

        #Control the top group to have the pressure target
        dy_top, Fv = Control_Top_PID(dict_algorithm, dict_sollicitations['Vertical_Confinement_Force'], dict_sample['L_g'])

        #Shear the top group and apply confinement force
        for grain in dict_sample['L_g'] :
            if grain.group == 'Top':
                grain.move_as_a_group(np.array([dict_sollicitations['Shear_velocity']*dict_algorithm['dt_DEM'], dy_top]), dict_algorithm['dt_DEM'])
        dict_sample['y_box_max'] = dict_sample['y_box_max'] + dy_top
        Shear_strain = Shear_strain + dict_sollicitations['Shear_velocity']*dict_algorithm['dt_DEM'] / Sample_height #Update shear strain

        #compute compacity
        Surface_g = 0
        for grain in dict_sample['L_g']:
            Surface_g = Surface_g + grain.surface

        #tracker
        dict_tracker['shear_L'].append(Shear_strain)
        dict_tracker['compacity_L'].append(Surface_g/((dict_sample['y_box_max']-dict_sample['y_box_min'])*(dict_sample['x_box_max']-dict_sample['x_box_min'])))
        dict_tracker['vertical_force_before_L'].append(sum_fy_top)
        dict_tracker['vertical_force_after_L'].append(Fv)
        dict_tracker['dy_top_L'].append(dy_top)
        dict_tracker['mu_sample_L'].append(mu_sample)

        #move solute out of grains
        if dict_algorithm['i_DEM'] % dict_algorithm['i_update_pf_solute'] == 0:
            #move pf for each grain with rbm
            simulation_report.tic_2nd_tempo()
            for i_grain in range(len(dict_sample['L_g'])) :
                #add bc
                dict_sample['L_g'][i_grain].move_grain_interpolation(dict_algorithm, dict_material, dict_sample)
                dict_sample['L_g'][i_grain].u_pf_interpolation = np.array([0,0])
                dict_sample['L_g'][i_grain].dtheta_pf_interpolation = 0
            #update etai
            for etai in dict_sample['L_etai']:
                etai.update_etai_M(dict_sample['L_g'])
            simulation_report.tac_2nd_tempo('Update phase fields')
            #move solute
            Owntools.Interpolate_solute_out_grains(dict_algorithm, dict_sample)

        #add pf simulation (new module) and pf->DEM and DEM->PF

        #debug print and plot
        if dict_algorithm['i_DEM'] % dict_algorithm['i_print_plot'] == 0:
            print('i_DEM',dict_algorithm['i_DEM'],': Confinement',int(100*Fv/dict_sollicitations['Vertical_Confinement_Force']),'% Shear',round(Shear_strain,4),'('+str(int(100*Shear_strain/dict_sollicitations['Shear_strain_target']))+' %)')
            if dict_algorithm['Debug_DEM'] :
                Owntools.Plot.Plot_Config_Sheared(dict_sample,dict_algorithm['i_DEM']) #change function here

        #Check stop conditions for DEM
        if Shear_strain >= dict_sollicitations['Shear_strain_target'] : #strain target reached
             DEM_loop_statut = False
        if dict_sample['L_g'] == []: #no more grains
            DEM_loop_statut = False
        if dict_algorithm['i_DEM'] >= dict_sollicitations['i_DEM_stop'] : #too many iterations
            DEM_loop_statut = False

    #plot total displacement field
    Owntools.Plot.Plot_total_U(dict_sample)
    #plot trackers
    Owntools.Plot.Plot_strain_compacity(dict_tracker)
    Owntools.Plot.Plot_strain_confinement(dict_tracker, dict_sollicitations)
    Owntools.Plot.Plot_strain_mu_sample(dict_tracker)
    Owntools.Plot.Plot_strain_dy_top(dict_algorithm, dict_tracker)

#-------------------------------------------------------------------------------

def Control_Top_PID(dict_algorithm, Force_target, L_g):
    """
    Control the upper wall to apply force.

    A PID corrector is applied.
        Input :
            an algorithm dictionnary (a dict)
            a confinement value (a float)
            a list of grain (a list)
        Output :
            the displacement of the top group (a float)
            a force applied on the top group before control (a float)
    """
    #compute vertical force applied on top group
    F = 0
    for grain in L_g :
        if grain.group == 'Top':
            F = F + grain.fy
    #compare with the target value
    error = F - Force_target #to have dy_top < 0 is F < Force_target
    #corrector
    ki = 0
    kd = 0
    dy_top = error * dict_algorithm['kp']
    #compare with maximum value
    if abs(dy_top) > dict_algorithm['dy_top_max'] :
        dy_top = np.sign(dy_top)*dict_algorithm['dy_top_max']

    return dy_top, F
