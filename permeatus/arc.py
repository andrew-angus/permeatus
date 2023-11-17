#!/bin/python3
# Author: Andrew Angus

import numpy as np
import csv
import matplotlib.pyplot as plt
import os
from importlib.resources import files
from permeatus.planar import *
import permeatus

# ABAQUS planar permeation class object
class arc(planar):

  # Initialisation arguments
  def __init__(self,extent=None,*args,**kwargs):

    # Call parent init
    super().__init__(*args,**kwargs)

    # Assign attributes
    self.extent = extent

    # Initialise derivative attributes

    # Set dispersed phase coefficients to zeros if none

  # Submit ABAQUS job
  def submit_job(self):

    # Update script template with variable inputs
    template = files(permeatus).joinpath("arc_nlayer.txt")
    with template.open() as i:
      with open("abaqus_script.py","w") as o:

        # Loop through template rows 
        for row in i:

          # Change desired lines
          if  row.startswith("dr = "):
            o.write(f'dr = {list(self.L)}\n')
          elif row.startswith("inner_r = "):
            o.write(f'inner_r = {self.r}\n')
          elif row.startswith("arc_extent = "):
            o.write(f'arc_extent = {self.extent}\n')
          elif row.startswith("D = "):
            o.write(f'D = {list(self.D)}\n')
          elif row.startswith("S = "):
            o.write(f'S = {list(self.S)}\n')
          elif row.startswith("N = "):
            o.write(f'N = {list(self.N)}\n')
          elif row.startswith("touts = "):
            o.write(f'touts = {list(self.touts)}\n')
          elif row.startswith("C0 = "):
            o.write(f'C0 = {self.C0}\n')
          elif row.startswith("C1 = "):
            o.write(f'C1 = {self.C1}\n')
          elif row.startswith("tstep = "):
            o.write(f'tstep = {self.tstep}\n')
          elif row.startswith("ncpu = "):
            o.write(f'ncpu = {self.ncpu}\n')

          # Write other lines unchanged
          else:
            o.write(row)

    # Remove output file to prevent appending to existing file
    try:
      os.system('rm C.csv')
      os.system('rm J.csv')
    except:
      pass

    # Submit created script
    print('Running ABAQUS...')
    os.system('abaqus cae noGui=abaqus_script.py')
    print('DONE')

  # Read field data output from abaqus csv file
  #TODO convert cartesian field read to polar coordinates
  #def read_field(self,target='C',targetdir=None):


  # Plot 1D solution
  #TODO update plot for 1d path through arc radius
  def plot_1d(self,target='C',showplot=True,timemask=None,plotlabels=None):
    pass

  # Linear algebra steady state solution
  #TODO pretty much the same but plotting needs to consider inner radius
  #def steady_state(self,y='C',plot=False,showplot=True,\
  #    plotlabel='steady-state'):

# Avoid divisions by zero
def mdiv(diver):
  return np.where(np.abs(diver) > np.finfo(np.float64).tiny, \
      diver, np.finfo(np.float64).tiny)
