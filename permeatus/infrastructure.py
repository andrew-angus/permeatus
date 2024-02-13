#!/bin/python3
# Author: Andrew Angus

import numpy as np
import csv
import matplotlib.pyplot as plt
import os
from importlib.resources import files
from permeatus.utils import *
import permeatus
import subprocess

# Infrastructure class object
class infrastructure:

  # Initialisation arguments
  def __init__(self,layers,L=None,touts=None,D=None,S=None,P=None,\
               C0=None,C1=None,p0=None,p1=None,\
               N=None,tstep=None,ncpu=None,\
               Dc=None,Sc=None,Pc=None,Vd_frac=None,AR=None,\
               Dd=None,Sd=None,Pd=None,model=None,solver='abaqus',\
               jobname='job',directory='.',verbose=True):

    # Defaults
    if ncpu is None:
      ncpu = 1
    if tstep is None:
      tstep = 0.001
    if C1 is None:
      C1 = 0.0

    #TODO Check validity of arguments 

    # Assign attributes
    self.layers = layers
    self.L = L
    self.D = D
    self.S = S
    self.P = P
    self.C0 = C0
    self.C1 = C1
    self.p0 = p0
    self.p1 = p1
    self.N = N
    self.touts = touts
    self.tstep = tstep
    self.ncpu = ncpu
    self.totL = np.sum(L)
    self.Dc = Dc
    self.Sc = Sc
    self.Pc = Pc
    self.Dd = Dd
    self.Sd = Sd
    self.Pd = Pd
    self.Vd_frac = Vd_frac
    self.AR = AR
    self.model = model
    self.solver = solver
    self.jobname = jobname
    self.directory = directory
    self.verbose = verbose

    # Initialise derivative attributes
    if self.touts is not None:
      self.frames = len(self.touts)+1
    self.field = {}
    self.J = None
    self.p = None
    self.C = None
    self.x = None
    self.xc = None
    self.Pavg = None
    self.Davg = None
    self.Savg = None
    self.P_upper = None
    self.D_upper = None
    self.S_upper = None
    self.P_lower = None
    self.D_lower = None
    self.S_lower = None

    # Set dispersed phase coefficients to zeros if none
    if Pd is None and Dd is None:
      Pd = np.zeros(layers)
      Dd = np.zeros(layers)
      Sd = np.zeros(layers)
    # Else get one from the other
    else:
      if self.Dd is not None:
        self.Pd = self.Dd*self.Sd
      elif self.Pd is not None:
        self.Dd = self.Pd
        self.Sd = np.ones_like(self.Dd)

    # Nielsen model if applicable
    if self.Vd_frac is not None:

      # Calculate Pc from DcSc, or vice versa
      if self.Dc is not None:
        self.Pc = self.Dc*self.Sc
      elif self.Pc is not None:
        self.Dc = self.Pc
        self.Sc = np.ones_like(self.Dc)

      # Evaluate effective coefficients
      if model == 'Nielsen':
        tortuosity = 1 + 0.5*self.AR*self.Vd_frac
        self.D = self.Dc/mdiv(tortuosity)
        self.S = self.Sc*(1-Vd_frac)
        self.P = self.D*self.S
      elif model == None:
        self.P = self.Pc
        self.D = self.Dc
        self.S = self.Sc
    else:
      # Calculate P from DS, or vice versa
      if self.D is not None:
        self.P = self.D*self.S
      elif self.P is not None:
        self.D = self.P
        self.S = np.ones_like(self.D)

    # Calculate one BC from the other
    if self.p0 is not None:
      self.C0 = self.p0*self.S[0]
    elif self.C0 is not None:
      self.p0 = self.C0/mdiv(self.S[0])
    if self.p1 is not None:
      self.C1 = self.p1*self.S[-1]
    elif self.C1 is not None:
      self.p1 = self.C1/mdiv(self.S[-1])

  # Setup planar problem using gmsh
  def __abaqus_input_file__(self):

    # Initialise
    gmsh.initialize()


    # Add model and set options
    gmsh.model.add(self.jobname)
    gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)
    gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)

    # Construct layers
    dx = self.totL
    dy = self.totL
    ticker = 0
    for i in range(self.layers):
      gmsh.model.occ.addRectangle(0,ticker,0,dx,self.L[i],tag=i)
      ticker += self.L[i]

    # Fragment overlappers
    out, pc = gmsh.model.occ.fragment([(2, i) for i in range(self.layers)], \
        [(2, i) for i in range(self.layers)])

    # Identify physical groups for material assignment
    gmsh.model.occ.synchronize()
    for i in range(self.layers):
      gmsh.model.addPhysicalGroup(2, [i], name=f"material{i}")

    # Define structured mesh with seeded edges
    eps = 1e-3*np.min(np.append(self.L,dx))
    ticker = 0
    for i in range(self.layers+1):
      ents = gmsh.model.getEntitiesInBoundingBox(-eps,-eps+ticker,-eps,\
          eps+dx,eps+ticker,eps,1)
      gmsh.model.mesh.setTransfiniteCurve(ents[0][1], 2)
      if i != self.layers:
        ents = gmsh.model.getEntitiesInBoundingBox(-eps,-eps+ticker,-eps,\
          eps,eps+ticker+self.L[i],eps,1)
        gmsh.model.mesh.setTransfiniteCurve(ents[0][1], self.N[i])
        ents = gmsh.model.getEntitiesInBoundingBox(-eps+dx,-eps+ticker,-eps,\
          eps+dx,eps+ticker+self.L[i],eps,1)
        gmsh.model.mesh.setTransfiniteCurve(ents[0][1], self.N[i])
        gmsh.model.mesh.setTransfiniteSurface(i)
        ticker += self.L[i]

    # Generate mesh
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.recombine()

    # Acquire boundary node sets
    gmsh.model.occ.synchronize()
    bottomnodes, topnodes, leftnodes, rightnodes = \
        boundary_nodes_2d(gmsh.model,dx,dy)

    # Write output and finalise
    write_abaqus_diffusion(self.D,self.S,self.C0,self.C1,self.touts,self.tstep,\
        bottomnodes,topnodes,leftnodes,rightnodes,self.jobname,PBC=False)
    gmsh.fltk.run()
    gmsh.finalize()

  # Submit FEA job (only abaqus for now
  def submit_job(self):

    if self.solver == 'abaqus':

      # Create input file
      self.__abaqus_input_file__()

      # Remove output file to prevent appending to existing file
      try:
        os.system('rm C.csv')
        os.system('rm J.csv')
      except:
        pass

      # Submit input file
      print('Running ABAQUS...')
      subprocess.call(\
          ['abaqus','interactive',f'job={self.jobname}',f'cpus={self.ncpu}'])

      # Obtain and update post-processing script
      postscript = files(permeatus).joinpath("data/abaqus_postscript.txt")
      with postscript.open() as i:
        with open("abaqus_postscript.py","w") as o:

          # Loop through template rows 
          for row in i:

            # Change jobname
            if  row.startswith("jobname = "):
              o.write(f"jobname = '{self.jobname}'\n")

            # Write other lines unchanged
            else:
              o.write(row)
      subprocess.call(['abaqus','cae','noGui=abaqus_postscript.py'])
      print('DONE')

  # Read field data output from abaqus csv file
  #TODO vector/tensor field reading
  def read_field(self,target='C',targetdir=None):

    # Target check
    targets = ['C','J']
    if target not in targets:
      raise Exception(f'target must be one of {targets}')

    # Initialise field-specific information
    self.field[target] = [None for i in range(self.frames)]
    field = self.field[target]
    if target == 'C':
      fieldkeys = ['CONC']
      fieldsize = 1
    elif target == 'J':
      fieldsize = 2
      fieldkeys = ['MFL-MFL1','MFL-MFL2']

    # Read into dictionary with csv module
    inc = -1; incc = -1
    fname = f'{target}.csv'
    labels = True
    if targetdir is not None:
      fname = os.path.join(targetdir,fname)
    with open(fname, newline='') as f:
      reader = csv.DictReader(f)
      for row in reader:
        # Strip whitespace from keys
        row = {i.strip():j for i,j in zip(row.keys(),row.values())}

        # Only store keys of interest
        keys = ['Frame','X','Y']+fieldkeys
        row = {i:row[i] for i in keys}

        # Split frame into increment and time and extract values
        splits = row['Frame'].split()
        increment = int(splits[1].strip(':'))
        time = float(splits[-1])

        # Check for new increment and initialise data dict
        if increment != inc:
          incc += 1
          field[incc] = {'t':time,'x':np.array([float(row['X'])]),\
              'y':np.array([float(row['Y'])]), \
              'data':np.array([[float(row[fieldkey]) for fieldkey in fieldkeys]])}
          inc = increment

        # Else append data
        else:
          field[incc]['x'] = np.r_[field[incc]['x'],float(row['X'])]
          field[incc]['y'] = np.r_[field[incc]['y'],float(row['Y'])]
          field[incc]['data'] = np.r_[field[incc]['data'],\
              np.array([[float(row[fieldkey]) for fieldkey in fieldkeys]])]

    # Store number of frames
    #self.field['frames'] = incc+1

  # Plot 1D solution
  def plot_1d(self,target='C',showplot=True,timemask=None,plotlabels=None):

    if timemask is None:
      timemask = [True for i in range(self.frames)]

    # Loop through frames
    for frame in range(1,self.frames):
      if timemask[frame]:
        # Identify path along part centreline in y-direction
        pathargs = np.ravel(np.argwhere(self.field['C'][frame]['x'] < 1e-10))
        y = self.field['C'][frame]['y'][pathargs]
        C = self.field['C'][frame]['data'][pathargs,0]

        # Sort by x-coordinate
        ysort = np.argsort(y)
        y = y[ysort]
        C = C[ysort]

        if self.layers > 1:
          # Get interfacial points
          Ly = np.zeros(self.layers+1)
          Ly[1:] += np.cumsum(self.L)
          iargsall = {}
          p = np.zeros(len(C)-self.layers+1)
          yp = np.unique(y)
          iargsold = [0,0]
          for i,j in enumerate(Ly[1:-1]):
            #print(y,j)
            iargs = np.argwhere(np.isclose(y,j)).flatten()#,atol=1e-14,rtol=1e-14)).flatten()
            iargsall[i] = iargs

            # Get interfacial pressure and swap points if not matching
            Cint = C[iargs]
            pint = Cint/self.S[i:i+2]
            if not np.isclose(pint[0],pint[1],atol=1e-8,rtol=1e-5):
              y[iargs[0]], y[iargs[1]] = y[iargs[1]], y[iargs[0]]
              C[iargs[0]], C[iargs[1]] = C[iargs[1]], C[iargs[0]]

            # Calculate final pressure array
            p[iargsold[1]-i:iargs[1]-i] = C[iargsold[1]:iargs[1]]/mdiv(self.S[i])
            iargsold = iargs
          p[iargsold[1]-(i+1):] = C[iargsold[1]:]/mdiv(self.S[-1])
        else:
          yp = y
          p = C/self.S[0]


        #TODO Plot either concentration or pressure
        if plotlabels is None or plotlabels[frame] is None:
          label = f"{self.field['C'][frame]['t']:0.3f} s"
        else:
          label = plotlabels[frame]
        if target == 'C':
          plt.plot(y,C,label=label)
        else:
          plt.plot(yp,p,label=label)

    # Finalise plotting
    plt.legend()
    plt.xlabel(r'$z$ [$m$]')
    if target == 'C':
      plt.ylabel(r'$C$ [mol$m^{-3}$]')
    else:
      plt.ylabel(r'$p$ [Pa]')
    if showplot:
      plt.show()

  # Function which reads xy report from abaqus
  def read_xy(self,fname='abaqus.rpt'):

    # Initialise data dict and labels logical
    self.xy = {}
    labels = True

    # Open report file and go through lines
    with open(fname) as f:
      for line in f:
        # Filter lines which are not of interest (anything but data and labels)
        if '             ' in line and 'Legend' not in line:
          # Split resulting lines into list of relevant labels/data
          entries = line.split(' ')
          entries = [i for i in entries if not (i == '' or i == '\n')]

          # Check for NoValue and strip new line characters
          skip = [False for i in range(len(entries))]
          for i in range(len(entries)):
            entries[i] = entries[i].rstrip()
            if entries[i] == 'NoValue':
              entries[i] = 0.0
              #skip[i] = True

          # Store labels
          if labels:
            self.xy = {i:np.empty(0) for i in entries}
            labels = False
          # Append data
          else:
            for i,label in enumerate(self.xy):
              if not skip[i]:
                self.xy[label] = np.append(self.xy[label],float(entries[i]))

  # Linear algebra steady state solution
  def steady_state(self,y='C',plot=False,showplot=True,\
      plotlabel='steady-state'):

    # Get linear coefficients relating pressure to molar flux
    k = self.P/mdiv(self.L)

    # Treat low layer number cases seperately
    # Evaluate pressure, concentration, and molar flux at relevant grid points
    if self.layers == 1:
      p = np.array([self.p0,self.p1])
      J = k[0]*(self.p0-self.p1)
      C = np.array([self.C0,self.C1])
      x = np.linspace(0,self.totL,2)
      xc = np.linspace(0,self.totL,2)
    elif self.layers == 2:
      p = np.zeros(3)
      p[0] = self.p0
      p[-1] = self.p1
      p[1] = k[0]*p[0]/mdiv(k[0]+k[1])
      J = k[0]*(p[0]-p[1])
      C = np.zeros(4)
      C[:2] = p[:2]*self.S[0]
      C[2:] = p[1:]*self.S[1]
      x = np.zeros_like(p)
      x[-1] = self.totL
      x[1] = self.L[0]
      xc = np.linspace(0,self.totL,4)
      xc[1:-1] = x[1]

    # n-layer solve
    else:
      # Populate internal pressure problem arrays
      b = np.zeros(self.layers-1)
      b[0] = self.p0*k[0]
      b[-1] = self.p1*k[-1]
      # coefficient matrix
      a = -k[1:-1]*(np.eye(self.layers-1,k=1)+np.eye(self.layers-1,k=-1)) + \
          (k[:-1]+k[1:])*np.eye(self.layers-1)

      # Solve for internal pressures
      p_int = np.linalg.solve(a,b)

      # Calculate molar flux
      J = k[0]*(self.p0-p_int[0])

      # Construct full pressure array
      p = np.zeros(self.layers+1)
      p[0] = self.p0
      p[-1] = self.p1
      p[1:-1] = p_int

      # Get x points for pressures
      x = np.zeros_like(p)
      x[1:] += np.cumsum(self.L)

      # Calculate concentrations from pressures
      # x points
      xc = np.zeros(2*self.layers)
      xc[-1] = x[-1]
      xc[1:-2:2] = x[1:-1]
      xc[2:-1:2] = x[1:-1]
      # concentrations
      C = np.zeros_like(xc)
      C[:-1:2] = self.S*p[:-1]
      C[1::2] = self.S*p[1:]

    # Optionally plot
    if plot:
      if y == 'C':
        plt.plot(xc,C,'-x',label=plotlabel)
        plt.ylabel(r'$C$ [mol$m^{-3}$]')
      else:
        plt.plot(x,p,'-x',label=plotlabel)
        plt.ylabel(r'$p$ [Pa]')
      plt.xlabel(r'$x$ [$m$]')
      if showplot:
        plt.show()

    # Store and return results
    self.x = x
    self.p = p
    self.J = J
    self.xc = xc
    self.C = C
    self.dp = np.diff(p)/self.L
    self.dC = np.array([(C[i+1]-C[i])/self.L[i//2] for i in range(0,len(C)-1,2)])

    if y == 'C':
      return xc, C, J
    else:
      return x, p, J

  # Calculate average steady state molar flux and store result
  def get_molar_flux(self,method='abaqus'):

    # Either calculate average of abaqus data or use steady-state 1D FEA
    if method == 'abaqus':
      if self.field is None:
        self.read_field('J')
      self.J = np.mean(self.field['J'][-1]['data'][:,1])
    elif method == 'steady-state':
      scrap,scrap,self.J = self.steady_state()

  #TODO calculate mass diffusion
  def get_mass_diffusion(self):
    pass

  #TODO calculate flowrate
  def get_flowrate(self):
    pass

  # For composite/crystalline structures, acquire average property bounds
  def get_eff_bnds(self,method='Wiener',set_eff='False'):

    # Wiener bounds, or absolute bounds by geometric and arithmetic means
    Vc_frac = 1-self.Vd_frac
    if method == 'Wiener':
      self.P_upper = Vc_frac*self.Pc + self.Vd_frac*self.Pd
      self.D_upper = Vc_frac*self.Dc + self.Vd_frac*self.Dd
      self.S_upper = Vc_frac*self.Sc + self.Vd_frac*self.Sd
      self.P_lower = 1/mdiv(Vc_frac/mdiv(self.Pc)+self.Vd_frac/mdiv(self.Pd))
      self.D_lower = 1/mdiv(Vc_frac/mdiv(self.Dc)+self.Vd_frac/mdiv(self.Dd))
      self.S_lower = 1/mdiv(Vc_frac/mdiv(self.Sc)+self.Vd_frac/mdiv(self.Sd))

    # Hashin-Strikman bounds, tighter than Wiener assuming structural isotropy
    elif method == 'HS':
      self.P_upper = np.zeros(self.layers)
      self.D_upper = np.zeros(self.layers)
      self.S_upper = np.zeros(self.layers)
      self.P_lower = np.zeros(self.layers)
      self.D_lower = np.zeros(self.layers)
      self.S_lower = np.zeros(self.layers)

      # Calculate upper and lower bounds (although which is which is yet to be
      # be determined
      Pcplus = self.Pc + self.Vd_frac/mdiv(1/mdiv(self.Pd-self.Pc) \
          +Vc_frac/mdiv(3*self.Pc))
      Pdplus = self.Pd + Vc_frac/(1/mdiv(self.Pc-self.Pd) \
          +self.Vd_frac/mdiv(3*self.Pd))
      Dcplus = self.Dc + self.Vd_frac/mdiv(1/mdiv(self.Dd-self.Dc) \
          +Vc_frac/mdiv(3*self.Dc))
      Ddplus = self.Dd + Vc_frac/mdiv(1/mdiv(self.Dc-self.Dd) \
          +self.Vd_frac/mdiv(3*self.Dd))
      Scplus = self.Sc + self.Vd_frac/mdiv(1/mdiv(self.Sd-self.Sc) \
          +Vc_frac/mdiv(3*self.Sc))
      Sdplus = self.Sd + Vc_frac/mdiv(1/mdiv(self.Sc-self.Sd) \
          +self.Vd_frac/mdiv(3*self.Sd))

      # Looping over layers required since relative magnitude \
      # of each phases coefficients important
      for i in range(self.layers):
        if self.Pc[i] > self.Pd[i]:
          self.P_upper[i] = Pcplus[i]
          self.P_lower[i] = Pdplus[i]
        else:
          self.P_upper[i] = Pdplus[i]
          self.P_lower[i] = Pcplus[i]
        if self.Dc[i] > self.Dd[i]:
          self.D_upper[i] = Dcplus[i]
          self.D_lower[i] = Ddplus[i]
        else:
          self.D_upper[i] = Ddplus[i]
          self.D_lower[i] = Dcplus[i]
        if self.Sc[i] > self.Sd[i]:
          self.S_upper[i] = Scplus[i]
          self.S_lower[i] = Sdplus[i]
        else:
          self.S_upper[i] = Sdplus[i]
          self.S_lower[i] = Scplus[i]

    # If desired, set effective coefficients as upper, lower, or avg
    if set_eff == 'upper':
      self.P = self.P_upper
      self.D = self.D_upper
      self.S = self.S_upper
    elif set_eff == 'lower':
      self.P = self.P_lower
      self.D = self.D_lower
      self.S = self.S_lower
    elif set_eff == 'avg':
      self.P = (self.P_lower+self.P_upper)*0.5
      self.D = (self.D_lower+self.D_upper)*0.5
      self.S = (self.S_lower+self.S_upper)*0.5

  # Get average coefficients for multiple layered system
  def get_avg_coeffs(self):

    # Avg coefficients are function of molar flux, system size and BCs
    #self.Davg = self.J*self.totL/mdiv(self.C0-self.C1)
    self.Davg = -self.J/np.sum(np.array([self.dC[i]*self.L[i]/self.totL \
        for i in range(self.layers)]))
    self.Pavg = -self.J/np.sum(np.array([self.dp[i]*self.L[i]/self.totL \
        for i in range(self.layers)]))
    #self.Pavg = self.J*self.totL/mdiv(self.p0-self.p1)
    self.Savg = self.Pavg/mdiv(self.Davg)

# Avoid divisions by zero
def mdiv(diver):
  return np.where(np.abs(diver) > np.finfo(np.float64).tiny, \
      diver, np.finfo(np.float64).tiny)

