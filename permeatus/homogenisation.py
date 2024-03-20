#!/bin/python3
# Author: Andrew Angus

import numpy as np
import csv
import matplotlib.pyplot as plt
import os
from importlib.resources import files
from permeatus.utils import *
from permeatus.layered1D import *
import permeatus
import subprocess
from scipy.optimize import brentq

# Homogenisation class object
class homogenisation(layered1D):

  # Initialisation arguments
  def __init__(self,materials,touts=None,D=None,S=None,P=None,\
               C0=None,C1=None,p0=None,p1=None,\
               tstep=None,ncpu=1,\
               vFrac=None,AR=None,solver='abaqus',\
               jobname='job',directory='.',verbose=True):

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
    self.vFrac = vFrac
    self.AR = AR
    self.solver = solver
    self.jobname = jobname
    self.directory = directory
    self.verbose = verbose

    # Initialise derivative attributes
    if self.touts is not None:
      self.frames = len(self.touts)+1
    self.field = {}
    self.Pavg = None
    self.Davg = None
    self.Savg = None
    self.nodeSets = [None for i in range(self.materials)]
    self.interfaceNodes = np.empty(0)

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

  def ebbermanmeshrectangular(self,r,lc):
    gmsh.initialize()

    # Add model and set options
    gmsh.model.add("random")
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
    #vFrac=1
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
    
    gmsh.fltk.run()
    
    gmsh.finalize()

                                      
  # Ebberman Paper Mesh  
  def ebbermanmesh(self,r,lc):
        
    V2 = self.vFrac[1]    
        
    w = ((2*3.14*0.26*0.26)/(1.732*V2))**0.5
    #h = 1.732*w
    h = w
    boxsize = h
    eps = 1e-3

    gmsh.initialize()
    gmsh.model.add("RVE1")
    
    # Bounding box
    gmsh.model.occ.addRectangle(0,0,0,w,h,tag=1)

    # Translations
    translationbase = np.array([-boxsize,0.0,boxsize])
    translations = np.array([[i,j] for i in translationbase for j in translationbase])

    # Add circles
    c = np.array([[0.0,0.0],[w/2,h/2]])
    nc = len(c)

    # Add circles with periodic wrapping
    boxdimtag,boxtag = periodic_disks(nc,c,gmsh.model,boxsize,r,eps)

    # Identify physical groups for material assignment
    ents = gmsh.model.getEntities(2)
    ents.remove(boxdimtag)
    gmsh.model.addPhysicalGroup(2, [boxtag], name="material0")
    gmsh.model.addPhysicalGroup(2, [i[1] for i in ents], name="material1")

    # Enforce periodic mesh on x-bounds
    periodic_mesh(gmsh.model,boxsize,eps)
    
    # Generate mesh
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMin",lc)
    gmsh.option.setNumber("Mesh.MeshSizeMax",lc)
    gmsh.model.mesh.generate(2)
    
    # Acquire boundary node sets
    gmsh.model.occ.synchronize()
    bottomnodes, topnodes, leftnodes, rightnodes = \
        boundary_nodes_2d(gmsh.model,w,h)

    # Write output and finalise
    write_abaqus_diffusion(self.D,self.S,self.C0,self.C1,self.touts,self.tstep,\
        bottomnodes,topnodes,leftnodes,rightnodes,self.jobname,PBC=True)
    
    gmsh.fltk.run()
    
    gmsh.finalize()

  # Create microstructure mesh by random insertion or 
  # Lubachevsky-Stillinger algorithm
  def cross_section_mesh(self,nc,r,minSpaceFac=0.1,maxMeshFac=0.1,\
      algorithm='LS',showmesh=True,seed=None):

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
          if bound_proximity_check_2d(newc[0],r,eps,boxsize):
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
      algorithm=algorithm,showmesh=showmesh)
    else:
      # Write output and finalise
      write_abaqus_diffusion(self.D,self.S,self.C0,self.C1,self.touts,self.tstep,\
          bottomnodes,topnodes,leftnodes,rightnodes,self.jobname,PBC=True)
      if showmesh:
        gmsh.fltk.run()
      gmsh.finalize()

  # Create mesh of Reuss bound setup
  def reuss_mesh(self,Nx=2,Ny=80):

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

    # Save material node sets
    for i in range(self.materials):
      self.nodeSets[i], coords = gmsh.model.mesh.getNodesForPhysicalGroup(2, i+1)

    # Interfacial nodes
    self.interfaceNodes = np.empty(0)
    for i in range(self.materials):
      for j in range(i):
        if i != j:
          for k in self.nodeSets[i]:
            for l in self.nodeSets[j]:
              if k == l:
                self.interfaceNodes = np.append(self.interfaceNodes,k)
    self.interfaceNodes = self.interfaceNodes.astype(np.intc)

    # Acquire boundary node sets
    gmsh.model.occ.synchronize()
    bottomnodes, topnodes, leftnodes, rightnodes = \
        boundary_nodes_2d(gmsh.model,boxsize,boxsize)

    # Write output and finalise
    write_abaqus_diffusion(self.D,self.S,self.C0,self.C1,self.touts,self.tstep,\
        bottomnodes,topnodes,leftnodes,rightnodes,self.jobname,PBC=True)
    gmsh.fltk.run()
    gmsh.finalize()

  # Create mesh of Voigt bound setup
  def voigt_mesh(self,Nx=20,Ny=40):

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
    gmsh.fltk.run()
    gmsh.finalize()

# test

  # Analytical Reuss bound
  def reuss_bound(self):

    self.P_eff = 1/np.sum(self.vFrac/self.P)
    self.D_eff = 1/np.sum(self.vFrac/self.D)
    self.S_eff = self.P_eff/self.D_eff

  # Analytical Voigt bound
  def voigt_bound(self):

    self.P_eff = np.sum(self.vFrac*self.P)
    #self.D_eff = np.sum(self.vFrac*self.D)
    #self.S_eff = self.P_eff/self.D_eff
    self.S_eff = np.sum(self.vFrac*self.S)
    self.D_eff = self.P_eff/self.S_eff

  # Hashin-Strikman upper bound
  def HS_upper_bound(self):

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
    #self.D_eff = res[1]
    self.S_eff = res[2]
    #self.S_eff = self.P_eff/self.D_eff
    self.D_eff = self.P_eff/self.S_eff

  # Hashin-Strikman lower bound
  def HS_lower_bound(self):

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
    #self.S_eff = res[2]
    self.S_eff = self.P_eff/self.D_eff
    #self.D_eff = self.P_eff/self.S_eff

  # Get prediction of Nielsen model
  def nielsen(self):

    # Check only 2 materials specified
    if self.materials != 2:
      raise Exception('method only implemented for 2 material system')

    # Evaluate effective coefficients
    tortuosity = 1 + 0.5*self.AR*self.vFrac[1]
    self.D_eff = self.D[0]/mdiv(tortuosity)
    self.S_eff = self.vFrac[0]*self.S[0]
    self.P_eff = self.D_eff*self.S_eff

  # Get prediction of Maxwell-Eucken model
  def maxwell_eucken(self):

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
    #self.S_eff = self.S[0]* \
    #    (2*self.S[0]+self.S[1]-2*(self.S[0]-self.S[1])*self.vFrac[1])/ \
    #    (2*self.S[0]+self.S[1]+(self.S[0]-self.S[1])*self.vFrac[1])
    #self.D_eff = self.P_eff/self.S_eff

  # Get prediction of Bruggeman model
  def bruggeman(self):

    # Check only 2 materials specified
    if self.materials != 2:
      raise Exception('method only implemented for 2 material system')

    # Evaluate effective coefficients
    Pterm = (3*self.vFrac[1]-1)*self.P[1]+(3*self.vFrac[0]-1)*self.P[0]
    self.P_eff = 0.25*(Pterm+np.sqrt(Pterm**2+8*self.P[0]*self.P[1]))
    Dterm = (3*self.vFrac[1]-1)*self.D[1]+(3*self.vFrac[0]-1)*self.D[0]
    self.D_eff = 0.25*(Dterm+np.sqrt(Dterm**2+8*self.D[0]*self.D[1]))
    self.S_eff = self.P_eff/self.D_eff
    #Sterm = (3*self.vFrac[1]-1)*self.S[1]+(3*self.vFrac[0]-1)*self.S[0]
    #self.S_eff = 0.25*(Sterm+np.sqrt(Sterm**2+8*self.S[0]*self.S[1]))
    #self.D_eff = self.P_eff/self.S_eff

  # Get prediction of Chen model
  def chen(self):

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
