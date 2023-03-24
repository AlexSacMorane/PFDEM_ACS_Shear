# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains ??.nt functions used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import math
import numpy as np

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Grain_Tempo_Polygonal:
  """
  A temporary grain used to generated an initial condition.
  """

#-------------------------------------------------------------------------------

  def __init__(self, grain_sphere, n_border):
    """Defining the grain.

        Input :
            a grain tempo (a tempo grain)
            a discretization of the grain (an int)
        Output :
            a grain tempo polygonal (a tempo polygonal grain)
    """
    L_border = []
    L_border_x = []
    L_border_y = []
    L_r = []
    L_theta_r = []
    #Build the border (for print)
    for i in range(n_border):
        theta = 2*math.pi*i/n_border
        p = np.array(grain_sphere.center) + np.array([grain_sphere.radius*math.cos(theta),grain_sphere.radius*math.sin(theta)])
        L_border.append(p)
        L_border_x.append(p[0])
        L_border_y.append(p[1])
        L_r.append(grain_sphere.radius)
        L_theta_r.append(theta)
    L_border.append(L_border[0])
    L_border_x.append(L_border_x[0])
    L_border_y.append(L_border_y[0])
    #save
    self.group = 'Current'
    self.image = None
    self.radius = grain_sphere.radius
    self.r_min = grain_sphere.radius
    self.r_max = grain_sphere.radius
    self.theta = 0
    self.rho_surf = grain_sphere.rho_surf
    self.surface = math.pi*grain_sphere.radius**2
    self.mass = math.pi*grain_sphere.radius**2*grain_sphere.rho_surf
    self.inertia = self.mass*grain_sphere.radius**2
    self.id = grain_sphere.id
    self.center = np.array(grain_sphere.center)
    self.l_border = L_border
    self.l_border_x = L_border_x
    self.l_border_y = L_border_y
    self.l_r = L_r
    self.l_theta_r = L_theta_r
    self.y = grain_sphere.y
    self.nu = grain_sphere.nu
    self.g = grain_sphere.g #shear modulus
    self.fx = 0
    self.fy = 0
    self.v = np.array([0,0])
    self.w = 0
    self.track_u = False

#-------------------------------------------------------------------------------

  def add_F(self, F, p_application):
      """
      Add a force to the grain.

        Input :
            itself (a grain_tempo)
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
            itself (a grain_tempo)
            a gravity value (a float)
        Output :
            Nothing, but attributes concerning the force applied are initialized (three floats)
      """
      self.fx = 0
      self.fy = -g*self.mass
      self.mz = 0

#-------------------------------------------------------------------------------

  def euler_semi_implicite(self,dt_DEM,factor):
    """
    Move the grain following a semi implicit euler scheme.

        Input :
            itself (a grain_tempo)
            a time step (a float)
            a factor to limite the displacement (a float)
        Output :
            Nothing, but the grain is moved
    """
    #translation
    a_i = np.array([self.fx,self.fy])/self.mass
    self.v = self.v + a_i*dt_DEM
    if np.linalg.norm(self.v) > self.radius*factor/dt_DEM: #limitation of the speed
        self.v = self.v * self.radius*factor/dt_DEM/np.linalg.norm(self.v)
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
            itself (a grain_tempo)
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
            itself (a grain_tempo)
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

#-------------------------------------------------------------------------------

class Grain_Image_Polygonal(Grain_Tempo_Polygonal):
  """
  An image grain used to generated an initial condition.
  """

#-------------------------------------------------------------------------------

  def __init__(self, real_grain, position):
    """Defining the image of a real grain.

        Input :
            itself (a grain_image)
            a real grain (a grain_tempo)
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
            a real grain (a grain_tempo)
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
