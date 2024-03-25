#!/bin/python3
# Author: Andrew Angus

"""Module for homogenisation of permeation in inhomogeneous systems.

The main component of this module is the homogenisation class, which inherits,
and builds on, functionality from the :mod:`permeatus.layered1D` class. 
The parent class contains routines for running ABAQUS simulations and 
post-processing into effective permeation coefficients. This class then 
contains methods to construct various finite element meshes representing 
inhomogenous permeation problems. It also contains analytical homogenisation 
models (effective medium theories). 

"""

import numpy as np
import gmsh
from permeatus.utils import *
from permeatus.layered1D import *
from scipy.optimize import brentq
from typing import Optional, Union
from typeguard import typechecked

# Exported objects
__all__ = ['homogenisation']

# Custom types
ArrayLike = Union[list,np.ndarray]

@typechecked
class homogenisation(layered1D):
  """Class for extracting effective properties from inhomogeneous systems.

  Attributes
  ----------
  materials
    Number of materials in the system.
  vFrac
    Volume fraction of each material
  D
    Diffusion coefficients for each material, 
    with suggested units: [mm\ :sup:`2`/hr].
  S
    Solubility coefficients for each material, 
    with suggested units: [nmol/mm\ :sup:`3`.MPa].
  P
    Permeability coefficients for each material, 
    with suggested units: [nmol/mm.hr.MPa].
  C0
    Concentration source boundary condition at bottom boundary,
    with suggested units: [nmol/mm\ :sup:`3`].
  C1
    Concentration sink boundary condition at top boundary, 
    with suggested units: [nmol/mm\ :sup:`3`].
  p0
    Pressure source boundary condition at bottom boundary, 
    with suggested units: [MPa].
  p1
    Pressure sink boundary condition at top boundary, 
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
  AR
    Aspect ratios for each material, if using a model in which aspect 
    ratio is accounted for.
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
  vFrac: np.ndarray
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
  AR: np.ndarray

  # Derived attributes
  field: dict
  frames: int
  P_eff: float
  D_eff: float
  S_eff: float

  def __init__(self, materials: int, vFrac: ArrayLike, 
          D: Optional[ArrayLike] = None, S: Optional[ArrayLike] = None, \
          P: Optional[ArrayLike] = None, C0: Optional[float] = None, 
          C1: Optional[float] = None, p0: Optional[float] = None, \
          p1: Optional[float] = None, touts: Optional[ArrayLike] = None, \
          tstep: Optional[float] = None, ncpu: int = 1, 
          jobname: str = 'job', verbose: bool =True, \
          AR: Optional[ArrayLike] = None):

    """Initialise class with arguments.

    Most arguments are optional, and only enforced if using functionality which
    requires them. If concentration problem parameters are specified, pressure 
    problem parameters will be calculated automatically from them, and vice-
    versa.

    """

    # Numpy array conversion
    if D is not None: D = np.array(D)  
    if S is not None: S = np.array(S)  
    if P is not None: P = np.array(P)  
    if touts is not None: touts = np.array(touts)  
    if vFrac is not None: touts = np.array(vFrac)  
    if AR is not None: touts = np.array(AR)  

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
    if len(D) != materials:
      raise Exception('D must be array like and same length as materials')
    if len(P) != materials:
      raise Exception('P must be array like and same length as materials')
    if len(S) != materials:
      raise Exception('S must be array like and same length as materials')
    if len(vFrac) != materials:
      raise Exception('vFrac must be array like and same length as materials')
    if AR is not None and len(AR) != materials:
      raise Exception('AR must be array like and same length as materials')
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
    self.vFrac = vFrac
    self.AR = AR

    # Initialise derivative attributes
    if self.touts is not None:
      self.frames = len(self.touts)+1
    else:
      self.frames = None
    self.field = {}
    self.P_eff = None
    self.D_eff = None
    self.S_eff = None



  def ebberman_mesh(self,r: float, lc: float, showMesh: bool = True):
    """Create mesh from Ebermann et. al paper

    Create microstructure RVE mesh given in :cite:t:`ebermannAnalytical2022`.

    Parameters
    ----------
    r
      Particle radius, suggested units: [mm]
    lc
      Mesh size control.
    showMesh
      Control whether to launch Gmsh GUI and show created mesh.

    """

    # Initialise
    gmsh.initialize()

    # Add model and set options
    gmsh.model.add(self.jobname)
    gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)

    # Bounding box
    boxwidth = 0.629e-3
    boxheight = 1.089e-3
    gmsh.model.occ.addRectangle(0,0,0,boxwidth,boxheight,tag=1)
    eps = 1e-5

    # Translations
    translationwidth = np.array([-boxwidth,0.0,boxwidth])
    translationheight = np.array([-boxheight,0.0,boxheight])
    translations = np.array([[i,j] for i in translationwidth for j in translationheight])

    # Add circles
    c = np.array([[0.0,0.0],[boxwidth/2,boxheight/2]])
    nc = len(c)

    # Add circles with periodic wrapping
    boxdimtag,boxtag = periodic_disks(nc,c,gmsh.model,boxwidth, boxheight, r,eps)
    
    # Identify physical groups for material assignment
    ents = gmsh.model.getEntities(2)
    ents.remove(boxdimtag)
    gmsh.model.addPhysicalGroup(2, [boxtag], name="material0")
    gmsh.model.addPhysicalGroup(2, [i[1] for i in ents], name="material1")
    
    # Enforce periodic mesh on x-bounds
    periodic_mesh(gmsh.model,boxwidth,eps)
    
    # Generate mesh
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMin",lc)
    gmsh.option.setNumber("Mesh.MeshSizeMax",lc)
    gmsh.model.mesh.generate(2)
    
    # Acquire boundary node sets
    gmsh.model.occ.synchronize()
    bottomnodes, topnodes, leftnodes, rightnodes = \
        boundary_nodes_2d(gmsh.model,boxwidth,boxheight)
    
    # Write output and finalise
    write_abaqus_diffusion(self.D,self.S,self.C0,self.C1,self.touts,self.tstep,\
        bottomnodes,topnodes,leftnodes,rightnodes,self.jobname,PBC=True)
    
    if showMesh:
      gmsh.fltk.run()
    
    gmsh.finalize()

                                      

  def cross_section_mesh(self, nc: int, r: float, minSpaceFac: float = 0.1, \
      maxMeshFac: float = 0.4, algorithm: str = 'LS', showMesh: bool = True, \
      seed: Optional[int] = None):

    """Create 2D fibre-reinforced composite perpendicular cross-section mesh

    Create fibrous composite microstructure mesh as 2D perpendicular cross-section.
    Available algorithms for microstructure creation are random insertion or
    Lubachevsky-Stillinger :cite:p:`lubachevskyGeometric1990`.

    Parameters
    ----------
    nc
      Number of circles, representing fibre cross-sections.
    r
      Fibre radius, suggested units: [mm]
    minSpaceFac
      Minimum spacing between fibres, as a factor of the fibre radius, to avoid
      meshing issues.
    maxMeshFac
      Maximum mesh size, as a factor of the fibre radius.
    algorithm
      Choose which algorithm to use, must be one of 'random' representing random
      insertion with acceptance-rejection, or 'LS' representing Lubachevsky-Stillinger.
    showMesh
      Control whether to launch Gmsh GUI and show created mesh.
    seed
      Integer seed for random number generator, which can allow reproduction of the same
      random mesh for testing.

    """

    # Check only 2 materials specified
    if self.materials != 2:
      raise Exception('cross_section_mesh only implemented for 2 material system')

    valids = ['LS','random']
    if algorithm not in valids:
      raise Exception(f'algorithm must be one of {valids}')

    # Random seed
    if seed is None:
      seed = np.random.randint(2e9)
    np.random.seed(seed)
		
    # Initialise
    gmsh.initialize()

    # Add model and set options
    gmsh.model.add(self.jobname)
    gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)
    gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)
    if not self.verbose:
      gmsh.option.setNumber("General.Terminal",0)

    # Microstructure specifications
    vfrac = self.vFrac[1] # Volume fraction

    # Bounding box - scale to desired volume fraction
    area0 = np.pi*r**2
    area1 = nc*area0
    area2 = area1/vfrac
    boxsize = np.sqrt(area2)
    gmsh.model.occ.addRectangle(0,0,0,boxsize,boxsize,tag=1)

    # Minimum spacing
    eps = r*minSpaceFac

    # Define buffers
    cbuff = 2*r + eps

    # Translations
    translationbase = np.array([-boxsize,0.0,boxsize])
    translations = np.array([[i,j] for i in translationbase for j in translationbase])

    ### Create microstructure ###

    # Lubachevsky-Stillinger algorithm
    if algorithm == 'LS':

      # LS algorithm
      # Randomly insert N point particles
      c = np.random.rand(nc,2)*np.array([[boxsize,boxsize]])

      # Randomly assign velocities
      v = np.random.rand(nc,2)*np.array([[2*boxsize,2*boxsize]])
      v[:,0] -= boxsize
      v[:,1] -= boxsize

      # Define algorithm variables
      t = 0 # Time
      rt = 0 # r(t)
      rf = r + eps/2 # Final radius, including minimum spacing
      h = rf # Radius growth rate
      tf = rf/h # Final time
      parthits = np.ones((nc,nc))*np.inf # Matrix of particle collision times
      niter = 0 # Iteration counter

      # Periodic wrapping of coordinates
      def wrap(point):
        point[0] -= np.floor(point[0]/boxsize)*boxsize
        point[1] -= np.floor(point[1]/boxsize)*boxsize
        return point

      # Collision time for particles with each other
      def ppColl(ci,cj,vi,vj,h,rt):
        cij = ci-cj
        vij = vi-vj
        A = np.dot(vij,vij) - 4*h**2
        B = np.dot(cij,vij) - 4*h*rt
        C = np.dot(cij,cij) - 4*rt**2
        quad = B**2-A*C
        if (B <= 0 or A < 0) and quad >=0:
          res =  1/A*(-B-np.sqrt(quad))
          return res
          if res > 0:
            return res
          else:
            return np.inf
        else:
          return np.inf

      # Simulate growth and collisions till circles have grown to r + eps/2
      advancing = True
      while advancing:

        # Get particle-particle collision times
        for i in range(nc):
          for j in range(i):

            # Get collision times with each periodic translation and take minimum
            parthits[i,j] = np.inf
            for translator in translations:
              tcoll = ppColl(c[i]+translator,c[j],v[i],v[j],h,rt)
              parthits[i,j] = np.minimum(parthits[i,j],tcoll)

        # Find overall minimum collision time
        colliders = np.unravel_index(parthits.argmin(), parthits.shape)
        time = parthits[colliders]
        
        # Check for final time
        if t + time > tf:

          # Final position update
          time = tf - t 
          for i in range(nc):
            c[i] = wrap(c[i]+v[i]*time)
          advancing = False

        else:
          # Propagate centers of colliding particles
          c[colliders[0]] = wrap(c[colliders[0]] + v[colliders[0]]*time)
          c[colliders[1]] = wrap(c[colliders[1]] + v[colliders[1]]*time)

          # Resolve collision
          dc = c[colliders[0]] - c[colliders[1]]
          dc[0] -= round(dc[0] / boxsize) * boxsize
          dc[1] -= round(dc[1] / boxsize) * boxsize
          u = dc/np.linalg.norm(dc)
          vpi = np.dot(u,v[colliders[0]])*u
          vti = v[colliders[0]] - vpi
          vpj = np.dot(u,v[colliders[1]])*u
          vtj = v[colliders[1]] - vpj
          v[colliders[0]] = vpj + vti
          v[colliders[1]] = vpi + vtj
          v[colliders[0]] += 2*u*h
          v[colliders[1]] -= 2*u*h
              
          # Propagate non-colliding centers
          for i in range(nc):
            if i != colliders[0] and i != colliders[1]:
              c[i] = wrap(c[i]+v[i]*time)
                
        # Update trackers
        t += time
        rt += h*time
        niter += 1

        # Check for issue
        if niter > 10000:
          print(niter,t,rt)
          raise Exception("Structure generation took too many iterations")

    # Random insertion algorithm
    else:

      # Loop over number of desired circles
      c = np.empty((0,2))
      for i in range(nc):

        # Continually try and insert circle within contstraints
        niter = 0
        reject = True
        while reject:

          # Randomly draw circle center
          newc = np.random.rand(1,2)*np.array([[boxsize,boxsize]])

          # For first circle just check minimum distance from edges and corners
          if bound_proximity_check_2d(newc[0],r,eps,boxsize,boxsize):
            if i == 0:
              reject = False

            # Afterwards check against distance between previous circles
            # and their periodic translations
            else:
              reject = False
              for center in c:
                mindist = np.inf
                for translator in translations:
                  dist = np.linalg.norm(newc[0]+translator-center)
                  mindist = np.minimum(dist,mindist)
                if mindist < cbuff:
                  reject = True
          niter += 1
          if niter > 5000:
              raise Exception("Structure generation took too many iterations")

        # Add accepted circle
        c = np.r_[c,newc]

    # Add circles with periodic wrapping
    boxdimtag,boxtag = periodic_disks(nc,c,gmsh.model,boxsize,boxsize,r,eps)
    print(boxdimtag,boxtag)

    # Identify physical groups for material assignment
    ents = gmsh.model.getEntities(2)
    ents.remove(boxdimtag)
    gmsh.model.addPhysicalGroup(2, [boxtag], name="material0")
    gmsh.model.addPhysicalGroup(2, [i[1] for i in ents], name="material1")

    # Enforce periodic mesh on x-bounds
    periodic_mesh(gmsh.model,boxsize,eps)

    # Generate mesh
    gmsh.option.setNumber("Mesh.MeshSizeMin",eps)
    gmsh.option.setNumber("Mesh.MeshSizeMax",r*maxMeshFac)
    gmsh.model.mesh.generate(2)

    # Acquire boundary node sets
    gmsh.model.occ.synchronize()
    bottomnodes, topnodes, leftnodes, rightnodes = \
        boundary_nodes_2d(gmsh.model,boxsize,boxsize)

    # Check for failed construction
    if len(leftnodes) != len(rightnodes):
      print("Warning: left node list not equal length to right nodes; rerunning with different seed")
      self.cross_section_mesh(nc=nc,r=r,minSpaceFac=minSpaceFac,maxMeshFac=maxMeshFac,\
      algorithm=algorithm,showMesh=showMesh)
    else:
      # Write output and finalise
      write_abaqus_diffusion(self.D,self.S,self.C0,self.C1,self.touts,self.tstep,\
          bottomnodes,topnodes,leftnodes,rightnodes,self.jobname,PBC=True)
      if showMesh:
        gmsh.fltk.run()
      gmsh.finalize()



  def reuss_mesh(self, Nx: int = 2, Ny: int = 80, showMesh: bool = True):
    """Create mesh whose analytical solution is the Reuss bound.

    Create a mesh of parallel material layers in the direction of flux.
    For detail see :cite:t:`auriaultHomogenization2010`.

    Parameters
    ----------
    Nx
      Number of cells in the x-direction.
    Ny
      Number of cells in the y-direction.
    showMesh
      Control whether to launch Gmsh GUI and show created mesh.

    """

    # Check only 2 materials specified
    if self.materials != 2:
      raise Exception('mesh only implemented for 2 material system')

    # Initialise
    gmsh.initialize()

    # Add model and set options
    gmsh.model.add(self.jobname)
    gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)
    gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)
    if not self.verbose:
      gmsh.option.setNumber("General.Terminal",0)

    # Layers
    vfrac = self.vFrac[1] # Volume fraction
    vfracfrac = Fraction(vfrac).limit_denominator()
    nlayer = vfracfrac.denominator
    ndlayer = vfracfrac.numerator
    layerid = np.random.choice(nlayer,size=ndlayer,replace=False)

    # Construct layers
    boxsize = 1
    dx = boxsize
    dy = boxsize/nlayer
    for i in range(nlayer):
        gmsh.model.occ.addRectangle(0,i*dy,0,dx,dy,tag=i)

    # Fragment overlappers
    out, pc = gmsh.model.occ.fragment([(2, i) for i in range(nlayer)], \
        [(2, i) for i in range(nlayer)])

    # Identify physical groups for material assignment
    gmsh.model.occ.synchronize()
    ents = gmsh.model.getEntities(2)
    dlayers = [(2,i) for i in layerid]
    ents = [i for i in ents if i not in dlayers]
    gmsh.model.addPhysicalGroup(2, [i[1] for i in ents], tag=1, name="material0")
    gmsh.model.addPhysicalGroup(2, [i[1] for i in dlayers], tag=2, name="material1")
    gmsh.model.occ.synchronize()

    # Define structured mesh with seeded edges
    eps = 1e-3/nlayer
    for i in range(2):
        ents = gmsh.model.getEntitiesInBoundingBox(-eps+i*dx,-eps,-eps,eps+i*dx,eps+boxsize,eps,1)
        for j in ents:
            gmsh.model.mesh.setTransfiniteCurve(j[1], Ny // nlayer + 1)
    for i in range(nlayer+1):
        ents = gmsh.model.getEntitiesInBoundingBox(-eps,-eps+i*dy,-eps,eps+dx,eps+i*dy,eps,1)
        for j in ents:
            gmsh.model.mesh.setTransfiniteCurve(j[1], Nx + 1)
        if i != nlayer:
            gmsh.model.mesh.setTransfiniteSurface(i)

    # Generate mesh
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.recombine()

    # Acquire boundary node sets
    gmsh.model.occ.synchronize()
    bottomnodes, topnodes, leftnodes, rightnodes = \
        boundary_nodes_2d(gmsh.model,boxsize,boxsize)

    # Write output and finalise
    write_abaqus_diffusion(self.D,self.S,self.C0,self.C1,self.touts,self.tstep,\
        bottomnodes,topnodes,leftnodes,rightnodes,self.jobname,PBC=True)
    if showMesh:
      gmsh.fltk.run()
    gmsh.finalize()



  def voigt_mesh(self, Nx: int = 20, Ny: int = 40, showMesh: bool = True):
    """Create mesh whose analytical solution is the Voigt bound.

    Create a mesh of material layers in series in the direction of flux.
    For detail see :cite:t:`auriaultHomogenization2010`.

    Parameters
    ----------
    Nx
      Number of cells in the x-direction.
    Ny
      Number of cells in the y-direction.
    showMesh
      Control whether to launch Gmsh GUI and show created mesh.

    """

    # Check only 2 materials specified
    if self.materials != 2:
      raise Exception('mesh only implemented for 2 material system')

    # Initialise
    gmsh.initialize()

    # Add model and set options
    gmsh.model.add(self.jobname)
    gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)
    gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)
    if not self.verbose:
      gmsh.option.setNumber("General.Terminal",0)

    # Layers
    vfrac = self.vFrac[1] # Volume fraction
    vfracfrac = Fraction(vfrac).limit_denominator()
    nlayer = vfracfrac.denominator
    ndlayer = vfracfrac.numerator
    layerid = np.random.choice(nlayer,size=ndlayer,replace=False)

    # Construct layers
    boxsize = 1
    dx = boxsize/nlayer
    dy = boxsize
    for i in range(nlayer):
      gmsh.model.occ.addRectangle(i*dx,0,0,dx,dy,tag=i)

    # Fragment overlappers
    out, pc = gmsh.model.occ.fragment([(2, i) for i in range(nlayer)], \
        [(2, i) for i in range(nlayer)])

    # Identify physical groups for material assignment
    gmsh.model.occ.synchronize()
    ents = gmsh.model.getEntities(2)
    dlayers = [(2,i) for i in layerid]
    ents = [i for i in ents if i not in dlayers]
    gmsh.model.addPhysicalGroup(2, [i[1] for i in ents], name="material0")
    gmsh.model.addPhysicalGroup(2, [i[1] for i in dlayers], name="material1")
    gmsh.model.occ.synchronize()

    # Define structured mesh with seeded edges
    eps = 1e-3/nlayer
    for i in range(2):
      ents = gmsh.model.getEntitiesInBoundingBox(-eps,-eps+i*dy,-eps,boxsize+eps,eps+i*dy,eps,1)
      for j in ents:
        gmsh.model.mesh.setTransfiniteCurve(j[1], Nx // nlayer + 1)
    for i in range(nlayer+1):
      ents = gmsh.model.getEntitiesInBoundingBox(-eps+i*dx,-eps,-eps,eps+i*dx,eps+dy,eps,1)
      for j in ents:
        gmsh.model.mesh.setTransfiniteCurve(j[1], Ny + 1)
      if i != nlayer:
        gmsh.model.mesh.setTransfiniteSurface(i)

    # Generate mesh
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.recombine()

    # Acquire boundary node sets
    gmsh.model.occ.synchronize()
    bottomnodes, topnodes, leftnodes, rightnodes = \
        boundary_nodes_2d(gmsh.model,boxsize,boxsize)

    # Write output and finalise
    write_abaqus_diffusion(self.D,self.S,self.C0,self.C1,self.touts,self.tstep,\
        bottomnodes,topnodes,leftnodes,rightnodes,self.jobname,PBC=True)
    if showMesh:
      gmsh.fltk.run()
    gmsh.finalize()



  def get_eff_coeffs(self):
    """Get effective coefficients of system by numerical averaging

    Get effective coefficients by numerical averaging. Simulation output at
    integration points is averaged by weighting by the volume of the integration
    point. The results are stored in the class attributes P_eff, D_eff, S_eff.

    """
    
    super().get_eff_coeffs(method='numerical')



  def reuss_bound(self):
    """Calculate analytical Reuss bound

    For detail see :cite:t`auriaultHomogenization2010`.

    """

    self.P_eff = 1/np.sum(self.vFrac/self.P)
    self.D_eff = 1/np.sum(self.vFrac/self.D)
    self.S_eff = self.P_eff/self.D_eff



  def voigt_bound(self):
    """Calculate analytical Voigt bound

    For detail see :cite:t:`auriaultHomogenization2010`.

    """

    self.P_eff = np.sum(self.vFrac*self.P)
    self.S_eff = np.sum(self.vFrac*self.S)
    self.D_eff = self.P_eff/self.S_eff



  def HS_upper_bound(self):
    """Calculate analytical Hashin-Strikman upper bound

    For detail see :cite:t:`auriaultHomogenization2010`.

    """

    # Check only 2 materials specified
    if self.materials != 2:
      raise Exception('method only implemented for 2 material system')

    # Determine high and low coefficients wrt permeation
    if self.P[1] > self.P[0]:
      hix = 1
      lox = 0
      vF = self.vFrac[0]
    else:
      hix = 0
      lox = 1
      vF = self.vFrac[1]

    # Calculate bound results
    res = []
    for i in [self.P,self.D,self.S]:
      hi, lo = i[hix], i[lox]
      res.append(hi + vF/(1/(lo-hi)+(1-vF)/(3*hi)))

    # Assign results
    self.P_eff = res[0]
    self.S_eff = res[2]
    self.D_eff = self.P_eff/self.S_eff



  def HS_lower_bound(self):
    """Calculate analytical Hashin-Strikman lower bound

    For detail see :cite:t:`auriaultHomogenization2010`.

    """

    # Check only 2 materials specified
    if self.materials != 2:
      raise Exception('method only implemented for 2 material system')

    # Determine high and low coefficients wrt permeation
    if self.P[1] > self.P[0]:
      hix = 1
      lox = 0
      vF = self.vFrac[0]
    else:
      hix = 0
      lox = 1
      vF = self.vFrac[1]

    # Calculate bound results
    res = []
    for i in [self.P,self.D,self.S]:
      hi, lo = i[hix], i[lox]
      res.append(lo + (1-vF)/(1/(hi-lo)+vF/(3*lo)))

    # Assign results
    self.P_eff = res[0]
    self.D_eff = res[1]
    self.S_eff = self.P_eff/self.D_eff

  def nielsen(self):
    """ Calculate homogenised coefficients by the Nielsen model.

    This model requires setting of the AR class attribute, and assumes
    an impermeable dispersed phase. See :cite:t:`prasadModeling2021` for detail.

    """

    # Check only 2 materials specified
    if self.materials != 2:
      raise Exception('method only implemented for 2 material system')

    # Evaluate effective coefficients
    tortuosity = 1 + 0.5*self.AR*self.vFrac[1]
    self.D_eff = self.D[0]/tortuosity
    self.S_eff = self.vFrac[0]*self.S[0]
    self.P_eff = self.D_eff*self.S_eff

  # Get prediction of Maxwell-Eucken model
  def maxwell_eucken(self):
    """ Calculate homogenised coefficients by the Maxwell-Eucken model.

    See :cite:t:`weiPredicting2018a` for detail.

    """

    # Check only 2 materials specified
    if self.materials != 2:
      raise Exception('method only implemented for 2 material system')

    # Evaluate effective coefficients
    self.P_eff = self.P[0]* \
        (2*self.P[0]+self.P[1]-2*(self.P[0]-self.P[1])*self.vFrac[1])/ \
        (2*self.P[0]+self.P[1]+(self.P[0]-self.P[1])*self.vFrac[1])
    self.D_eff = self.D[0]* \
        (2*self.D[0]+self.D[1]-2*(self.D[0]-self.D[1])*self.vFrac[1])/ \
        (2*self.D[0]+self.D[1]+(self.D[0]-self.D[1])*self.vFrac[1])
    self.S_eff = self.P_eff/self.D_eff

  # Get prediction of Bruggeman model
  def bruggeman(self):
    """ Calculate homogenised coefficients by the Bruggeman model.

    See :cite:t:`weiPredicting2018a` for detail.

    """

    # Check only 2 materials specified
    if self.materials != 2:
      raise Exception('method only implemented for 2 material system')

    # Evaluate effective coefficients
    Pterm = (3*self.vFrac[1]-1)*self.P[1]+(3*self.vFrac[0]-1)*self.P[0]
    self.P_eff = 0.25*(Pterm+np.sqrt(Pterm**2+8*self.P[0]*self.P[1]))
    Dterm = (3*self.vFrac[1]-1)*self.D[1]+(3*self.vFrac[0]-1)*self.D[0]
    self.D_eff = 0.25*(Dterm+np.sqrt(Dterm**2+8*self.D[0]*self.D[1]))
    self.S_eff = self.P_eff/self.D_eff

  # Get prediction of Chen model
  def chen(self):
    """ Calculate homogenised coefficients by the Chen model.

    See :cite:t:`chenHomogenization2002` for detail.

    """

    # Check only 2 materials specified
    if self.materials != 2:
      raise Exception('method only implemented for 2 material system')

    # Evaluate effective coefficients
    def Proot(P_eff):
      f = P_eff/self.P[0]
      x = self.P[1]/self.P[0]
      return self.vFrac[0]**2*((1-x)/(f-x))**2*f - 1
    self.P_eff = brentq(Proot,self.P[0],self.P[1])
    def Droot(D_eff):
      f = D_eff/self.D[0]
      x = self.D[1]/self.D[0]
      return self.vFrac[0]**2*((1-x)/(f-x))**2*f - 1
    self.D_eff = brentq(Droot,self.D[0],self.D[1])
    self.S_eff = self.P_eff/self.D_eff

  # Shuts off redundant inherited steady state method
  def steady_state(self):
    """Obsolete inherited method; not implemented.
    """
    raise NotImplementedError

  # Shuts off redundant inherited plot_1d method
  def plot_1d(self):
    """Obsolete inherited method; not implemented.
    """
    raise NotImplementedError
