# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains functions used in the simulation to define grains.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import math
import numpy as np

#Own
import Owntools

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Grain:
  """
  A polygonal used to simulate the grain.
  """

#-------------------------------------------------------------------------------

  def __init__(self, grain_tempo):
    """
    Defining the grain.

        Input :
            a grain tempo (a tempo grain)
        Output :
            a real grain (a grain)
    """
    self.group = grain_tempo.group
    self.image = None
    self.radius = grain_tempo.radius
    self.r_min = grain_tempo.r_min
    self.r_max = grain_tempo.r_max
    self.theta = grain_tempo.theta
    self.rho_surf = grain_tempo.rho_surf
    self.surface = grain_tempo.surface
    self.mass = grain_tempo.mass
    self.inertia = grain_tempo.inertia
    self.id = grain_tempo.id
    self.center = grain_tempo.center.copy()
    self.l_border = grain_tempo.l_border.copy()
    self.l_border_x = grain_tempo.l_border_x.copy()
    self.l_border_y = grain_tempo.l_border_y.copy()
    self.l_r = grain_tempo.l_r.copy()
    self.l_theta_r = grain_tempo.l_theta_r.copy()
    self.y = grain_tempo.y
    self.nu = grain_tempo.nu
    self.g = grain_tempo.g #shear modulus
    self.fx = 0
    self.fy = 0
    self.v = np.array([0,0])
    self.mz = 0
    self.w = 0
    self.track_u = False

#-------------------------------------------------------------------------------

  def add_F(self, F, p_application):
      """
      Add a force to the grain.

        Input :
            itself (a grain)
            a force applied (a 1 x 2 numpy array)
            a application point (a 1 x 2 numpy array)
        Output :
            Nothing, but attributes are updated (three floats)
      """
      self.fx = self.fx + F[0]
      self.fy = self.fy + F[1]
      v1 = np.array([p_application[0]-self.center[0], p_application[1]-self.center[1], 0])
      v2 = np.array([F[0], F[1], 0])
      self.mz = self.mz + np.cross(v1,v2)[2]

#-------------------------------------------------------------------------------

  def init_F_control(self,g):
      """
      Initialize the force applied to the grain.

      A gravity is assumed.

        Input :
            itself (a grain)
            a gravity value (a float)
        Output :
            Nothing, but attributes concerning the force applied are initialized (three floats)
      """
      self.fx = 0
      self.fy = -g*self.mass
      self.mz = 0

#-------------------------------------------------------------------------------

  def euler_semi_implicite(self,dt_DEM):
    """
    Move the grain following a semi implicit euler scheme.

        Input :
            itself (a grain)
            a time step (a float)
        Output :
            Nothing, but the grain is moved
    """
    #translation
    a_i = np.array([self.fx,self.fy])/self.mass
    self.v = self.v + a_i*dt_DEM
    for i in range(len(self.l_border)):
        self.l_border[i] = self.l_border[i] + self.v*dt_DEM
        self.l_border_x[i] = self.l_border_x[i] + self.v[0]*dt_DEM
        self.l_border_y[i] = self.l_border_y[i] + self.v[1]*dt_DEM
    self.center = self.center + self.v*dt_DEM
    if self.track_u :
        self.total_ux = self.total_ux + self.v[0]*dt_DEM
        self.total_uy = self.total_uy + self.v[1]*dt_DEM

    #rotation
    dw_i = self.mz/self.inertia
    self.w = self.w + dw_i*dt_DEM
    self.theta = self.theta + self.w*dt_DEM

    for i_theta_r in range(len(self.l_theta_r)) :
        theta_r = self.l_theta_r[i_theta_r]
        theta_r = theta_r + self.w*dt_DEM
        while theta_r >= 2*math.pi:
            theta_r = theta_r - 2*math.pi
        while theta_r < 0 :
            theta_r = theta_r + 2*math.pi
        self.l_theta_r[i_theta_r] = theta_r

    for i in range(len(self.l_border)):
        p = self.l_border[i] - self.center
        Rot_Matrix = np.array([[math.cos(self.w*dt_DEM), -math.sin(self.w*dt_DEM)],
                               [math.sin(self.w*dt_DEM),  math.cos(self.w*dt_DEM)]])
        p = np.dot(Rot_Matrix,p)
        self.l_border[i] = p + self.center
        self.l_border_x[i] = p[0] + self.center[0]
        self.l_border_y[i] = p[1] + self.center[1]

#-------------------------------------------------------------------------------

  def move_as_a_group(self, U, dt_DEM):
    """
    Move the grain in a group defined.

        Input :
            itself (a grain)
            a displacement (a 2 x 1 numpy array)
        Output :
            Nothing, but the grain is moved
    """
    #translation
    self.v = np.array([U[0]/dt_DEM, U[1]/dt_DEM])
    for i in range(len(self.l_border)):
        self.l_border[i] = self.l_border[i] + U
        self.l_border_x[i] = self.l_border_x[i] + U[0]
        self.l_border_y[i] = self.l_border_y[i] + U[1]
    self.center = self.center + U
    if self.track_u :
        self.total_ux = self.total_ux + U[0]
        self.total_uy = self.total_uy + U[1]

    #rotation
    self.w = 0

#-------------------------------------------------------------------------------

  def is_group(self, ymin, ymax, name_group):
      """
      Check if a grain is in a determined group by comparing center y-coordinate with two limits.

        Input :
            itself (a grain)
            two y limits (two floats)
            a name (a string)
        Output :
            a Boolean and the attribut group is updated (a string)

      """
      if ymin <= self.center[1] and self.center[1] <= ymax:
          self.group = name_group
          #stop the grain
          self.v = np.array([0,0])
          self.w = 0
          return True
      else :
          return False

#---------------------------------------------------------------------------

  def build_etai_M(self,dict_material,dict_sample):
        '''
        Build the phase field for one grain.

        A cosine profile is assumed (see https://mooseframework.inl.gov/source/ics/SmoothCircleIC.html).

            Input :
                itself (a grain)
                a material dictionnary (a dictionnary)
                a sample dictionnary (a dictionnary)
            Output :
                Nothing but the grain gets a new attribute (a n_y x n_x numpy array)
        '''
        #initilization
        self.etai_M = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))))

        #extract a spatial zone
        x_min = min(self.l_border_x)-dict_material['w']
        x_max = max(self.l_border_x)+dict_material['w']
        y_min = min(self.l_border_y)-dict_material['w']
        y_max = max(self.l_border_y)+dict_material['w']

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

        for l in range(i_y_min,i_y_max+1):
            for c in range(i_x_min,i_x_max+1):
                y = dict_sample['y_L'][l]
                x = dict_sample['x_L'][c]
                p = np.array([x,y])
                r = np.linalg.norm(self.center - p)
                #look for the radius on this direction
                if p[1]>self.center[1]:
                    theta = math.acos((p[0]-self.center[0])/np.linalg.norm(self.center-p))
                else :
                    theta= 2*math.pi - math.acos((p[0]-self.center[0])/np.linalg.norm(self.center-p))
                L_theta_R_i = list(abs(np.array(self.l_theta_r)-theta))
                R = self.l_r[L_theta_R_i.index(min(L_theta_R_i))]
                #build etai_M
                self.etai_M[-1-l][c] = Owntools.Cosine_Profile(R,r,dict_material['w'])

#-------------------------------------------------------------------------------

class Grain_Image(Grain):
  """
  An image grain used to onsider the periodic boundary conditions.
  """

#-------------------------------------------------------------------------------

  def __init__(self, real_grain, position):
    """
    Defining the image of a real grain.

        Input :
            itself (a grain_image)
            a real grain (a grain)
        Output :
            Nothing, but an image grain is generated (a grain_image)
    """
    real_grain.image = self
    self.real = real_grain
    self.position = position
    self.group = real_grain.group
    self.radius = real_grain.radius
    self.r_min = real_grain.radius
    self.r_max = real_grain.radius
    self.theta = real_grain.theta
    self.rho_surf = real_grain.rho_surf
    self.surface = real_grain.surface
    self.mass = real_grain.mass
    self.inertia = real_grain.inertia
    self.id = real_grain.id
    self.center = np.array(real_grain.center.copy())
    self.l_border = real_grain.l_border.copy()
    self.l_border_x = real_grain.l_border_x.copy()
    self.l_border_y = real_grain.l_border_y.copy()
    self.l_r = real_grain.l_r.copy()
    self.l_theta_r = real_grain.l_theta_r.copy()
    self.y = real_grain.y
    self.nu = real_grain.nu
    self.g = real_grain.g
    self.fx = real_grain.fx
    self.fy = real_grain.fy
    self.v = real_grain.v.copy()
    self.w = real_grain.w

#-------------------------------------------------------------------------------

  def translation(self, U):
    """Translate an image grain depending on the tempo grain.

        Input :
            itself (a grain_image)
            a real grain (a grain)
        Output :
            Nothing, but an image grain is generated (a grain_image)
    """
    self.center = self.real.center.copy() + U
    self.l_border = []
    self.l_border_x = []
    self.l_border_y = []
    for i in range(len(self.real.l_border)):
        self.l_border.append(self.real.l_border[i].copy() + U)
        self.l_border_x.append(self.real.l_border[i][0] + U[0])
        self.l_border_y.append(self.real.l_border[i][1] + U[1])
    #update data
    self.v = self.real.v.copy()
    self.w = self.real.w
    self.theta = self.real.theta
    self.l_r = self.real.l_r.copy()
    self.l_theta_r = self.real.l_theta_r.copy()
