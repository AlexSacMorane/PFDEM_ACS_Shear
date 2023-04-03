# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define a new class.
The new class is about the unconserved variables used in phase-field simulation.
"""

#-------------------------------------------------------------------------------
#Libs
#-------------------------------------------------------------------------------

import numpy as np
import random

#Own
import Report
import Grain

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Etai:

    #-------------------------------------------------------------------------------
    def __init__(self, ID, L_ig, L_g):
        '''
        Defining the phase variables etai.

        One eta represents several grains.

            Input :
                an id (an integer)
                a list of grain id represented by the variable (a list)
                the list of the grains (a list)
            Output :
                an eta (an eta)
        '''
        self.id = ID
        self.l_ig = L_ig
        self.update_etai_M(L_g)
        for i_grain in self.l_ig :
            L_g[i_grain].etai = self.id

    #-------------------------------------------------------------------------------

    def update_etai_M(self,L_g):
        '''
        Sum variables field of all the grain assigned to this eta.

            Input :
                itself (an eta)
                the list of grains (a list)
            Output :
                Nothing, but the phase field attribute is updated (a nx x ny numpy array)
        '''
        self.etai_M = np.array(L_g[self.l_ig[0]].etai_M.copy())
        for i in range(1,len(self.l_ig)):
            self.etai_M = self.etai_M + np.array(L_g[self.l_ig[i]].etai_M.copy())

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def etai_distribution(dict_algorithm, dict_sample, simulation_report, group):
    '''
    Assign grains to etai if there are in group defined.

    A minimal distance between grains with same eta is computed from the factor variable

        Input :
            an algorithm dictionnary (a dict)
            a sample dictionnary (a dict)
            a simulation report (a report)
            a group name (a string)
        Output :
            Nothing, but the sample dictionnary is updated with the list of the etai (a list)
    '''
    #xtraction of grains in the group
    L_g_extract = []
    for grain in dict_sample['L_g']:
        if grain.group == group:
            L_g_extract.append(grain)
    #initialisation
    for grain in L_g_extract:
        grain.ig_near_L = []
        grain.eta_near_L = []
    #look for grains near
    for i_grain1 in range(1,len(L_g_extract)):
        grain1 = L_g_extract[i_grain1]
        for i_grain2 in range(0,i_grain1):
            grain2 = L_g_extract[i_grain2]
            if np.linalg.norm(grain1.center - grain2.center) < dict_algorithm['factor_distribution_etai'] * (grain1.r_max+grain2.r_max):
                grain1.ig_near_L.append(i_grain2)
                grain2.ig_near_L.append(i_grain1)
    #first try
    n_etai_L = [1]
    L_etai = [len(dict_sample['L_etai'])]
    L_ig_etai = [[0]]
    L_g_extract[0].id_eta = len(dict_sample['L_etai'])
    for grain in L_g_extract:
        if len(dict_sample['L_etai']) in grain.ig_near_L:
            grain.eta_near_L.append(len(dict_sample['L_etai']))
    for i_grain in range(1,len(L_g_extract)):
        grain = L_g_extract[i_grain]
        etai_defined = False
        etai = len(dict_sample['L_etai'])
        while etai < len(L_etai) and not etai_defined:
            if etai not in grain.eta_near_L:
                grain.id_eta = etai
                etai_defined = True
                n_etai_L[etai] = n_etai_L[etai] + 1
                L_ig_etai[etai].append(i_grain)
                for grain2 in L_g_extract:
                    if i_grain in grain2.ig_near_L:
                        grain2.eta_near_L.append(etai)
            etai = etai + 1
        if etai == len(L_etai) and not etai_defined:
            grain.id_eta = etai
            L_etai.append(etai)
            n_etai_L.append(1)
            L_ig_etai.append([i_grain])
            for grain2 in L_g_extract:
                if i_grain in grain2.ig_near_L:
                    grain2.eta_near_L.append(etai)
    #adaptation (try to have the same number of grain assigned to all eta)
    n_etai_mean = np.mean(n_etai_L)
    #check the quality
    adaptation_done = True
    for n_etai in n_etai_L :
        if n_etai_mean-2 < n_etai and n_etai < n_etai_mean+2:
            adaptation_done = False
    adaptation_i = 0
    while not adaptation_done :
        adaptation_i = adaptation_i + 1
        L_ig_over =  L_ig_etai[n_etai_L.index(max(n_etai_L))]
        id_g_to_work = L_ig_over[random.randint(0,len(L_ig_over)-1)]
        i_etai_over = n_etai_L.index(max(n_etai_L))
        etai_over = L_etai[i_etai_over]
        i_etai_under = n_etai_L.index(min(n_etai_L))
        etai_under = L_etai[i_etai_under]
        grain = L_g_extract[id_g_to_work]
        if etai_under not in grain.eta_near_L:
            grain.id_eta = etai_under
            n_etai_L[i_etai_over] = n_etai_L[i_etai_over] - 1
            n_etai_L[i_etai_under] = n_etai_L[i_etai_under] + 1
            L_ig_etai[i_etai_over].remove(id_g_to_work)
            L_ig_etai[i_etai_under].append(id_g_to_work)
            for grain2 in L_g_extract:
                if id_g_to_work in grain2.ig_near_L:
                    grain2.eta_near_L.remove(etai_over)
                    grain2.eta_near_L.append(etai_under)
        #check the quality
        adaptation_done = True
        for n_etai in n_etai_L :
            if n_etai_mean-2 < n_etai and n_etai < n_etai_mean+2:
                adaptation_done = False
        if adaptation_i > len(L_g_extract):
            adaptation_done = True

    #Create etai
    L_etai = []
    for i in range(len(L_ig_etai)) :
        etai = Etai(len(dict_sample['L_etai'])+i, L_ig_etai[i], dict_sample['L_g'])
        L_etai.append(etai)

    #update the dict
    dict_sample['L_etai'] = dict_sample['L_etai'] + list(L_etai)
    simulation_report.write_and_print(f"{len(dict_sample['L_etai'])} phase variables used for group "+group+".\n",f"{len(dict_sample['L_etai'])} phase variables used for group "+group+".")

#-------------------------------------------------------------------------------
#Not used
#-------------------------------------------------------------------------------

def PFtoDEM_Multi(FileToRead,dict_algorithm,dict_material,dict_sample):
    '''
    Read file from MOOSE simulation to reconstruct the phase field of the grain.

        Input :
            the name of the file to read (a string)
            an algorithm dictionnary (a dictionnary)
            a material dictionnary (a dictionnary)
            a sample dictionnary (a dictionnary)
        Output :
            Nothing but the grain gets an updated attribute (a n_y x n_x numpy array)
    '''
    #--------------------------------------------------------------------------
    #Global parameters
    #---------------------------------------------------------------------------

    L_etai_M = np.array(np.zeros((len(dict_sample['L_etai']),len(dict_sample['y_L']),len(dict_sample['x_L']))))

    id_L = None
    eta_selector_len = len('        <DataArray type="Float64" Name="eta')
    end_len = len('        </DataArray>')
    XYZ_selector_len = len('        <DataArray type="Float64" Name="Points"')
    data_jump_len = len('          ')

    for i_proc in range(dict_algorithm['np_proc']):

        L_Work = [[], #X
                  []] #Y
        for etai in dict_sample['L_etai']:
            L_Work.append([]) #etai

    #---------------------------------------------------------------------------
    #Reading file
    #---------------------------------------------------------------------------

        f = open(f'{FileToRead}_{i_proc}.vtu','r')
        data = f.read()
        f.close
        lines = data.splitlines()

        #iterations on line
        for line in lines:

            if line[0:eta_selector_len] == '        <DataArray type="Float64" Name="eta':
                id_L = 1 + int(line[eta_selector_len])

            elif line[0:XYZ_selector_len] == '        <DataArray type="Float64" Name="Points"':
                id_L = 0

            elif (line[0:end_len] == '        </DataArray>' or  line[0:len('          <InformationKey')] == '          <InformationKey') and id_L != None:
                id_L = None

            elif id_L != None :
                if line[0:data_jump_len] == '          ' and id_L >= 2: #Read etai
                    line = line[data_jump_len:]
                    c_start = 0
                    for c_i in range(0,len(line)):
                        if line[c_i]==' ':
                            c_end = c_i
                            L_Work[id_L].append(float(line[c_start:c_end]))
                            c_start = c_i+1
                    L_Work[id_L].append(float(line[c_start:]))

                elif line[0:data_jump_len] == '          ' and id_L == 0: #Read [X, Y, Z]
                    line = line[data_jump_len:]
                    XYZ_temp = []
                    c_start = 0
                    for c_i in range(0,len(line)):
                        if line[c_i]==' ':
                            c_end = c_i
                            XYZ_temp.append(float(line[c_start:c_end]))
                            if len(XYZ_temp)==3:
                                L_Work[0].append(XYZ_temp[0])
                                L_Work[1].append(XYZ_temp[1])
                                XYZ_temp = []
                            c_start = c_i+1
                    XYZ_temp.append(float(line[c_start:]))
                    L_Work[0].append(XYZ_temp[0])
                    L_Work[1].append(XYZ_temp[1])

        #Adaptating data and update of etai_M
        for i in range(len(L_Work[0])):
            #Interpolation method
            L_dy = []
            for y_i in dict_sample['y_L'] :
                L_dy.append(abs(y_i - L_Work[1][i]))
            L_dx = []
            for x_i in dict_sample['x_L'] :
                L_dx.append(abs(x_i - L_Work[0][i]))
            for j in range(2,len(L_Work)):
                L_etai_M[j-2][-1-list(L_dy).index(min(L_dy))][list(L_dx).index(min(L_dx))] = L_Work[j][i]

    #---------------------------------------------------------------------------
    #Transmit data to grains
    #---------------------------------------------------------------------------

    for grain in dict_sample['L_g']:
        grain.ExtractPF_from_Eta(L_etai_M, dict_algorithm, dict_material, dict_sample)
        grain.geometric_study(dict_sample)
