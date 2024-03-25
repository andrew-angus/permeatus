#!/bin/python3
# Author: Andrew Angus

"""Module for modelling of permeation in 1D layered systems.

The main component of this module is the layered1D class, which contains 
methods and attributes for modelling 1D layered systems using finite element
analysis. Functionality includes an interface with ABAQUS, as well as an
alternative steady-state solver. It also includes post-processing, with
plotting of solutions and extracton of effective homogenised coefficients.

"""

# Module imports
import numpy as np
import csv
import matplotlib.pyplot as plt
import os
from importlib.resources import files
from permeatus.utils import *
import permeatus
import subprocess
from typing import Optional, Union, Tuple
from typeguard import typechecked

# Exported objects
__all__ = ['layered1D']

# Custom types
ArrayLike = Union[list,np.ndarray]

@typechecked
class layered1D:
  """Class for modelling 1D layered systems.

  This class contains attributes and methods for modelling 1D layered
  systems, with either ABAQUS for time-dependent solutions, or with
  a steady-state finite element solver. Capablities include modelling of
  both pressure and concentration driven problems, as well as plotting of 
  solutions, and extraction of effective properties.

  Attributes
  ----------
  materials
    Number of material layers in the system.
  L
    1D lengths of each layer in direction of permeation, 
    with suggested units: [mm].
  D
    Diffusion coefficients for each layer, 
    with suggested units: [mm\ :sup:`2`/hr].
  S
    Solubility coefficients for each layer, 
    with suggested units: [nmol/mm\ :sup:`3`.MPa].
  P
    Permeability coefficients for each layer, 
    with suggested units: [nmol/mm.hr.MPa].
  C0
    Concentration source boundary condition at first layer,
    with suggested units: [nmol/mm\ :sup:`3`].
  C1
    Concentration sink boundary condition at last layer, 
    with suggested units: [nmol/mm\ :sup:`3`].
  p0
    Pressure source boundary condition at first layer, 
    with suggested units: [MPa].
  p1
    Pressure sink boundary condition at last layer, 
    with suggested units: [MPa].
  touts
    Solution output times, if using ABAQUS,
    with suggested units: [hr].
  tstep
    Simulation timestep, if using ABAQUS,
    with suggested units: [hr].
  ncpu
    Number of CPUs to utilise, if using ABAQUS.
  jobname
    Name of job, if using ABAQUS.
  verbose
    Boolean flag which switches verbose output on or off.
  totL
    Total length of layered system.
  field
    Dictionary of solution data at integration points, which can contain 
    pressure, concentration, molar flux, and pressure/concentration 
    gradient data for each timeframe. Additionally, integration point 
    volume data can be stored.
  frames
    Integer number of frames, corresponding to number of touts.
  D_eff
    Effective diffusion coefficient, derived from solution.
  S_eff
    Effective solubility coefficient, derived from solution.
  P_eff
    Effective permeability coefficient, derived from solution.

  """

  # Init attributes
  materials: int
  L: np.ndarray
  D: np.ndarray
  S: np.ndarray
  P: np.ndarray
  C0: float
  C1: float
  p0: float
  p1: float
  touts: np.ndarray
  tstep: float
  ncpu: int
  jobname: str
  verbose: bool

  # Derived attributes
  totL: float
  field: dict
  frames: int
  P_eff: float
  D_eff: float
  S_eff: float

  def __init__(self, materials: int, L: ArrayLike, 
      D: Optional[ArrayLike] = None, S: Optional[ArrayLike] = None, \
      P: Optional[ArrayLike] = None, C0: Optional[float] = None, 
      C1: Optional[float] = 0.0, p0: Optional[float] = None, \
      p1: Optional[float] = 0.0, touts: Optional[ArrayLike] = None, \
      tstep: Optional[float] = None, ncpu: Optional[int] = 1, 
      jobname: Optional[str] = 'job', verbose: bool =True):

    """Initialise class with arguments.

    Most arguments are optional, and only enforced if using functionality which
    requires them. If concentration problem parameters are specified, pressure 
    problem parameters will be calculated automatically from them, and vice-
    versa.

    """
    # Numpy array conversion
    L = np.array(L)
    if D is not None: D = np.array(D)  
    if S is not None: S = np.array(S)  
    if P is not None: P = np.array(P)  
    if touts is not None: touts = np.array(touts)  

    # Calculate P from DS, or vice versa
    if D is not None and S is not None:
      P = D*S
    elif P is not None:
      if D is None and S is None:
        D = P
        S = np.ones_like(P)
      elif D is None:
        D = P/S
      else:
        S = P/D
    else:
      raise Exception('Either P or D & S must be specified')

    # Calculate one BC from the other
    if p0 is not None:
      C0 = p0*S[0]
    elif C0 is not None:
      p0 = C0/S[0]
    else:
      raise Exception('One of p0 or C0 must be specified.')
    if p1 is not None:
      C1 = p1*S[-1]
    elif C1 is not None:
      p1 = C1/S[-1]
    else:
      p1 = C1 = 0

    # Check validity of arguments 
    if materials < 1:
      raise Exception('materials argument must be an integer greater than 0')
    if len(L) != materials:
      raise Exception('L must be array like and same length as materials')
    if len(D) != materials:
      raise Exception('D must be array like and same length as materials')
    if len(P) != materials:
      raise Exception('P must be array like and same length as materials')
    if len(S) != materials:
      raise Exception('S must be array like and same length as materials')
    if tstep is not None and tstep <= 0.0:
      raise Exception('tstep must be float greater than 0')
    if ncpu < 1:
      raise Exception('ncpu argument must be an integer greater than 0')
    if p0 < 0.0:
      raise Exception('p0 must be >= 0')
    if p1 < 0.0:
      raise Exception('p1 must be >= 0')
    if C0 < 0.0:
      raise Exception('C0 must be >= 0')
    if C1 < 0.0:
      raise Exception('C1 must be >= 0')

    # Assign attributes
    self.materials = materials
    self.L = L
    self.D = D
    self.S = S
    self.P = P
    self.C0 = C0
    self.C1 = C1
    self.p0 = p0
    self.p1 = p1
    self.touts = touts
    self.tstep = tstep
    self.ncpu = ncpu
    self.jobname = jobname
    self.verbose = verbose

    # Initialise derivative attributes
    self.totL = np.sum(L)
    if self.touts is not None:
      self.frames = len(self.touts)+1
    else:
      self.frames = None
    self.field = {}
    self.P_eff = None
    self.D_eff = None
    self.S_eff = None



  def layered_mesh(self, N: ArrayLike):
    """Create mesh using Gmsh for solving layered system in Abaqus.

    Parameters
    ----------

    N
      Integer number of computational cells assigned to each layer, 
      if modelling with ABAQUS.

    """

    # Check for N, touts, tstep
    try:
      N = np.array(N)
    except:
      raise Exception('Argument N must be given as list of integers'+ \
          ' same length as number of materials')
    if type(N) != np.ndarray or type(N[0]) != int or len(N) != self.materials:
      raise Exception('Argument N must be given as list of integers'+ \
          ' same length as number of materials')
    if self.touts is None:
      raise Exception('Attribute touts must be specified.')
    if self.tstep is None:
      raise Exception('Attribute tstep must be specified.')

    # Initialise
    gmsh.initialize()

    # Add model and set options
    gmsh.model.add(self.jobname)
    gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)
    gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)
    if not self.verbose:
      gmsh.option.setNumber("General.Terminal",0)

    # Construct materials
    dx = self.totL
    dy = self.totL
    ticker = 0
    for i in range(self.materials):
      gmsh.model.occ.addRectangle(0,ticker,0,dx,self.L[i],tag=i)
      ticker += self.L[i]

    # Fragment overlappers
    out, pc = gmsh.model.occ.fragment([(2, i) for i in range(self.materials)], \
        [(2, i) for i in range(self.materials)])

    # Identify physical groups for material assignment
    gmsh.model.occ.synchronize()
    for i in range(self.materials):
      gmsh.model.addPhysicalGroup(2, [i], tag=i, name=f"material{i}")

    # Define structured mesh with seeded edges
    eps = 1e-3*np.min(np.append(self.L,dx))
    ticker = 0
    for i in range(self.materials+1):
      ents = gmsh.model.getEntitiesInBoundingBox(-eps,-eps+ticker,-eps,\
          eps+dx,eps+ticker,eps,1)
      gmsh.model.mesh.setTransfiniteCurve(ents[0][1], 2)
      if i != self.materials:
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



  def submit_job(self):
    """Submit Abaqus job.
    
    Mesh creation should be conducted prior to job submission.

    The output files C.csv, J.csv, and V.csv are produced, with concentration,
    flux, and integration point volumes respectively. These are produced from a
    python script submitted to Abaqus cae which processes the .odb output database.
    The script template is located in permeatus/data/abaqus_postscript.py, relative
    to the package root directory.

    """

    # Remove output file to prevent appending to existing file
    try:
      self.field = {}
      subprocess.call(['rm',"C.csv"],stderr=subprocess.DEVNULL)
      subprocess.call(['rm',"V.csv"],stderr=subprocess.DEVNULL)
      subprocess.call(['rm',"J.csv"],stderr=subprocess.DEVNULL)
    except:
      pass

    # Submit input file
    if self.verbose:
      print('Running ABAQUS...')
      subprocess.call(\
          ['abaqus','interactive',f'job={self.jobname}',f'cpus={self.ncpu}'])
    else:
      subprocess.call(\
          ['abaqus','interactive',f'job={self.jobname}',f'cpus={self.ncpu}'],\
          stdout=subprocess.DEVNULL, \
          stderr=subprocess.STDOUT)

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

    # Pass script to abaqus cae to produce csv output files
    if self.verbose:
      subprocess.call(['abaqus','cae','noGui=abaqus_postscript.py'])
      print('DONE')
    else:
      subprocess.call(['abaqus','cae','noGui=abaqus_postscript.py'],\
      stdout=subprocess.DEVNULL, \
      stderr=subprocess.STDOUT)



  def read_field(self,target: str):
    """Read field data output from Abaqus csv file

    Process Abaqus output csv file into class attribute field dictionary.
    Field dictionary has the following layout: 
    field[<field>][<frame>] = {'x': <x co-ordinates>, 'y': <y co-ordinates>,
    'data': <field data>, 'material': <material number>}
    
    The integration volume field is only stored for a single frame 0, as it
    will be unchanged through time.

    Must have run submit_job() to produce the required csv files.

    Parameters
    ----------

    target
      Specify target field; must be one of 'C', 'J', or 'V'.

    """

    # Target check
    targets = ['C','J','V']
    if target not in targets:
      raise Exception(f'target must be one of {targets}')

    # Initialise field-specific information
    if target != 'V':
      self.field[target] = [None for i in range(self.frames)]
    else:
      self.field[target] = [None]
    field = self.field[target]
    if target == 'C':
      fieldkeys = ['CONC']
      fieldsize = 1
    elif target == 'V':
      fieldkeys = ['IVOL']
      fieldsize = 1
    elif target == 'J':
      fieldsize = 2
      fieldkeys = ['MFL-MFL1','MFL-MFL2']

    # Read into dictionary with csv module
    inc = -1; incc = -1
    fname = f'{target}.csv'
    labels = True
    with open(fname, newline='') as f:
      reader = csv.DictReader(f)
      for row in reader:
        # Strip whitespace from keys
        row = {i.strip():j for i,j in zip(row.keys(),row.values())}

        # Only store keys of interest
        keys = ['Frame','Material Name','X','Y']+fieldkeys
        row = {i:row[i] for i in keys}

        # Split frame into increment and time and extract values
        splits = row['Frame'].split()
        increment = int(splits[1].strip(':'))
        time = float(splits[-1])

        # Get material if integration point output
        if target in ['C','J','V']:
          material = int(row['Material Name'].strip('MATERIAL'))
        else:
          material = -1

        # Check for new increment and initialise data dict
        if increment != inc:
          incc += 1
          field[incc] = {'t':time,'x':np.array([float(row['X'])]),\
              'y':np.array([float(row['Y'])]),\
              'material':np.array([material]), \
              'data':np.array([[float(row[fieldkey]) for fieldkey in fieldkeys]])}
          inc = increment

        # Else append data
        else:
          field[incc]['x'] = np.r_[field[incc]['x'],float(row['X'])]
          field[incc]['y'] = np.r_[field[incc]['y'],float(row['Y'])]
          field[incc]['material'] = np.r_[field[incc]['material'],material]
          field[incc]['data'] = np.r_[field[incc]['data'],\
              np.array([[float(row[fieldkey]) for fieldkey in fieldkeys]])]



  def get_eff_coeffs(self,method='numerical'):
    """Get effective coefficients of system by desired method

    Get effective coefficients by either numerical averaging or by 
    analytical solution to layered system (Reuss bound solution). The
    results are stored in the class attributes P_eff, D_eff, S_eff.

    Parameters
    ----------

    method
      Desired method; one of 'numerical' or 'analytical'

    """

    if method == 'numerical':

      # Check if concentration and pressure gradient data calculated
      if self.field == {} or 'C' not in self.field.keys() or \
          'J' not in self.field.keys() or 'p' not in self.field.keys() or \
          'grad' not in self.field['C'][-1].keys() or \
          'grad' not in self.field['p'][-1].keys():
        self.get_gradC()
        self.get_gradp()

      # Check if volume data available
      if 'V' not in self.field.keys():
        self.read_field('V')

      self.D_eff = -self.V_mean(self.field['J'][-1]['data'][:,1])/\
          self.V_mean(self.field['C'][-1]['grad'][:,1])
      self.P_eff = -self.V_mean(self.field['J'][-1]['data'][:,1])/ \
          self.V_mean(self.field['p'][-1]['grad'][:,1])
      self.S_eff = self.P_eff/self.D_eff

    elif method == 'analytical':

      self.P_eff = 1/np.sum(self.L/(self.totL*self.P))
      self.D_eff = 1/np.sum(self.L/(self.totL*self.D))
      self.S_eff = self.P_eff/self.D_eff



  def get_P_eff(self):
    """Get effective permeability of system by numerical averaging.
    
    The result is stored in the P_eff class attribute.

    """
    # Check if concentration and pressure gradient data calculated
    if self.field == {} or \
        'J' not in self.field.keys() or 'p' not in self.field.keys() or \
        'grad' not in self.field['p'][-1].keys():
      self.get_gradp()

    # Check if volume data available
    if 'V' not in self.field.keys():
      self.read_field('V')

    self.P_eff = -self.V_mean(self.field['J'][-1]['data'][:,1])/ \
        self.V_mean(self.field['p'][-1]['grad'][:,1])



  def V_mean(self,field: ArrayLike) -> float:
    """Return integration point volume weighted mean of field.
    
    Parameters
    ----------
    field
      The field to be averaged by integration point volume weighting.  

    Returns
    -------

    float
      Volume-weighted average of field



    """
    sumfieldV = np.sum(field*self.field['V'][-1]['data'][:,0])
    sumV = np.sum(self.field['V'][-1]['data'])
    return sumfieldV/sumV



  def get_gradC(self):
    """Calculate concentration gradient at integration points

    Calculate concentration gradient field from flux and concentration solutions,
    via Fick's first law. Store result in dictionary under 'grad' key.

    """

    # Read fields if not already read
    if 'J' not in self.field.keys():
      self.read_field('J')
    if 'C' not in self.field.keys():
      self.read_field('C')

    # Loop through frames and calculate gradient by Fick's first law
    for i in range(self.frames):
      D = np.array([self.D[j] for j in self.field['J'][i]['material']])[:,None]
      self.field['C'][i]['grad'] = -self.field['J'][i]['data']/D



  def get_gradp(self):
    """Calculate pressure gradient at integration points

    Calculate pressure gradient field from flux and pressure solutions,
    via Darcy's law. Store result in dictionary under 'grad' key.

    """

    # Get pressure field if non-existent
    if 'p' not in self.field.keys():
      self.get_p()
    if 'J' not in self.field.keys():
      self.read_field('J')

    # Loop through frames and calculate gradient by Darcy's first law
    for i in range(self.frames):
      P = np.array([self.P[j] for j in self.field['J'][i]['material']])[:,None]
      self.field['p'][i]['grad'] = -self.field['J'][i]['data']/P



  def get_p(self):
    """Calculate pressure field at integration points

    Calculate pressure field from concentration solution and solubilities.
    Store result in field attribute.

    """
    # Get concentration field if non-existent
    if 'C' not in self.field.keys():
      self.read_field('C')

    # Loop through frames and calculate pressure via solubility coefficient
    self.field['p'] = [None for i in range(self.frames)]
    for i in range(self.frames):
      self.field['p'][i] = {}
      self.field['p'][i]['x'] = self.field['C'][i]['x']
      self.field['p'][i]['y'] = self.field['C'][i]['y']
      self.field['p'][i]['material'] = self.field['C'][i]['material']
      S = np.array([self.S[j] for j in self.field['C'][i]['material']])
      self.field['p'][i]['data'] = -self.field['C'][i]['data']/S



  # Plot 1D solution
  def plot_1d(self, plotTarget: str = 'C', showPlot: bool = True, \
      timemask: Optional[ArrayLike] = None):
    """Plot 1D layered solution from numerical data

    Numerical fields for the target should be read previous to calling this method.

    Parameters
    ----------

    plotTarget
      One of 'C' or 'p', to determine whether concentration or pressure solutions 
      are plotted.
    showPlot
      Control whether to show plot immediately or delay (useful to compile multiple
      plots in one figure).
    timemask
      Boolean array the same length as the number of frames, which controls whether
      to plot that frame.

    """

    # Default time mask is true for all frames
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

        if self.materials > 1:
          # Get interfacial points
          Ly = np.zeros(self.materials+1)
          Ly[1:] += np.cumsum(self.L)
          iargsall = {}
          p = np.zeros(len(C)-self.materials+1)
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
            p[iargsold[1]-i:iargs[1]-i] = C[iargsold[1]:iargs[1]]/self.S[i]
            iargsold = iargs
          p[iargsold[1]-(i+1):] = C[iargsold[1]:]/self.S[-1]
        else:
          yp = y
          p = C/self.S[0]


        # Plot either concentration or pressure
        if plotlabels is None or plotlabels[frame] is None:
          label = f"{self.field['C'][frame]['t']:0.3f} s"
        if target == 'C':
          plt.plot(y,C,label=label)
        elif target == 'p':
          plt.plot(yp,p,label=label)

    # Finalise plotting
    plt.legend()
    plt.xlabel(r'x')
    if target == 'C':
      plt.ylabel(r'C')
    else:
      plt.ylabel(r'p')
    if showplot:
      plt.show()



  def steady_state(self, plot: bool = False, plotTarget: str = 'C', \
      showPlot: bool = True) \
      -> Tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray,float]:
    """1D finite element steady state solution

    Get steady-state solution with native 1D finite element solver. The method 
    calculates the pressure and concentration solutions at layer boundaries, as well
    as the scalar flux solution. The concentration solution will be dual-valued at 
    internal boundaries. Optional plotting is controlled by function arguments.

    Parameters
    ----------

    plot
      Control whether to plot solution.
    plotTarget
      One of 'C' or 'p', to determine whether concentration or pressure solutions 
      are plotted.
    showPlot
      Control whether to show plot immediately or delay (useful to compile multiple
      plots in one figure).

    Returns
    -------

    x
      Spatial points of pressure solution
    p
      Pressure solution
    xc
      Spatial points of concentration solution
    C
      Concentration solution
    J
      Scalar flux solution

    """

    # Get linear coefficients relating pressure to molar flux
    k = self.P/self.L

    # Treat low layer number cases seperately
    # Evaluate pressure, concentration, and molar flux at relevant grid points
    if self.materials == 1:
      p = np.array([self.p0,self.p1])
      J = k[0]*(self.p0-self.p1)
      C = np.array([self.C0,self.C1])
      x = np.linspace(0,self.totL,2)
      xc = np.linspace(0,self.totL,2)
    elif self.materials == 2:
      p = np.zeros(3)
      p[0] = self.p0
      p[-1] = self.p1
      p[1] = k[0]*p[0]/(k[0]+k[1])
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
      b = np.zeros(self.materials-1)
      b[0] = self.p0*k[0]
      b[-1] = self.p1*k[-1]
      # coefficient matrix
      a = -k[1:-1]*(np.eye(self.materials-1,k=1)+np.eye(self.materials-1,k=-1)) + \
          (k[:-1]+k[1:])*np.eye(self.materials-1)

      # Solve for internal pressures
      p_int = np.linalg.solve(a,b)

      # Calculate molar flux
      J = k[0]*(self.p0-p_int[0])

      # Construct full pressure array
      p = np.zeros(self.materials+1)
      p[0] = self.p0
      p[-1] = self.p1
      p[1:-1] = p_int

      # Get x points for pressures
      x = np.zeros_like(p)
      x[1:] += np.cumsum(self.L)

      # Calculate concentrations from pressures
      # x points
      xc = np.zeros(2*self.materials)
      xc[-1] = x[-1]
      xc[1:-2:2] = x[1:-1]
      xc[2:-1:2] = x[1:-1]
      # concentrations
      C = np.zeros_like(xc)
      C[:-1:2] = self.S*p[:-1]
      C[1::2] = self.S*p[1:]

    # Optionally plot
    if plot:
      if plotTarget == 'C':
        plt.plot(xc,C,'-x',label=plotlabel)
        plt.ylabel('C')
      else:
        plt.plot(x,p,'-x',label=plotlabel)
        plt.ylabel('p')
      plt.xlabel(r'x')
      if showplot:
        plt.show()

    # Return results
    return x, xc, p, C, J

