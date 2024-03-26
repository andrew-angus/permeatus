#!/bin/python3
# Author: Andrew Angus

"""Utility functions for creating meshes and setting up ABAQUS permeation simulations.

These utility functions mostly utilise gmsh objects, and introduce functionality to either
aid in mesh creation, or to set up ABAQUS permeation simulations.

"""

import numpy as np
import gmsh
from pandas import unique
import fileinput as fi
from typing import Optional, Tuple

# Functions which are exported
__all__ = ['boundary_nodes_2d','bound_proximity_check_2d','periodic_copy',\
    'periodic_disks','periodic_mesh','nodeset']



def boundary_nodes_2d(m: gmsh.model, dx: float, dy: float) \
    -> Tuple[np.ndarray[float], np.ndarray[float],\
    np.ndarray[float],np.ndarray[float]]:
  """Get boundary node sets in 2D box setups

  Given a Gmsh model, and bounding box x and y dimensions, extract
  the node sets for each boundary.

  Arguments
  ---------

  m
    Gmsh model; shortcut for gmsh.model.
  dx
    Bounding box x dimensions.
  dy
    Bounding box y dimensions.

  Returns
  -------

  bottomnodes
    Node set for bottom boundary.
  topnodes
    Node set for top boundary.
  leftnodes
    Node set for left boundary.
  rightnodes
    Node set for right boundary.


  """

  # Set up empty node and coordinate lists for each boundary
  leftnodes = np.empty(0)
  leftcoords = np.empty(0)
  rightnodes = np.empty(0)
  rightcoords = np.empty(0)
  topnodes = np.empty(0)
  topcoords = np.empty(0)
  bottomnodes = np.empty(0)
  bottomcoords = np.empty(0)

  # Get boundary entities
  boundary = m.getBoundary(m.getEntities(2))

  # Loop through boundary lines
  for dim, tag in boundary:

    # Get nodes and coordinates for this line
    nodeTags, coord, parametricCoord = \
        m.mesh.getNodes(dim, np.abs(tag), includeBoundary=True)
    coords = np.reshape(coord, (len(coord)//3, 3))

    # Determine which boundary the line resides on
    # and store nodes/coordinates in associated list
    main_coord = np.argmax(np.abs(np.sum(np.diff(coords,axis=0),axis=0)))
    if main_coord == 0:
      if np.isclose(coord[1],0.0):
        bottomnodes = np.r_[bottomnodes,nodeTags]
        bottomcoords = np.r_[bottomcoords, coords[:,main_coord]]
      elif np.isclose(coord[1],dy):
        topnodes = np.r_[topnodes,nodeTags]
        topcoords = np.r_[topcoords, coords[:,main_coord]]
    elif main_coord == 1:
      if np.isclose(coord[0],0.0):
        leftnodes = np.r_[leftnodes,nodeTags]
        leftcoords = np.r_[leftcoords, coords[:,main_coord]]
      elif np.isclose(coord[0],dx):
        rightnodes = np.r_[rightnodes,nodeTags]
        rightcoords = np.r_[rightcoords, coords[:,main_coord]]

  # Sort node lists by coordinates
  sorts = np.argsort(bottomcoords)
  bottomcoords = bottomcoords[sorts]
  bottomnodes = unique(bottomnodes[sorts]).astype(np.intc)
  sorts = np.argsort(topcoords)
  topcoords = topcoords[sorts]
  topnodes = unique(topnodes[sorts]).astype(np.intc)
  sorts = np.argsort(leftcoords)
  leftcoords = leftcoords[sorts]
  leftnodes = unique(leftnodes[sorts][1:-1]).astype(np.intc)
  sorts = np.argsort(rightcoords)
  rightcoords = rightcoords[sorts]
  rightnodes = unique(rightnodes[sorts][1:-1]).astype(np.intc)

  # Return node lists
  return bottomnodes, topnodes, leftnodes, rightnodes



def bound_proximity_check_2d(c: np.ndarray[float], r: float, eps: float, \
    dx: float, dy: float) -> bool:
  """Check for disk proximities to 2D bounding box boundaries.

  Check given disks outer edges are a minimum distance from bounding box
  boundaries.

  Arguments
  ---------

  c
    x and y coordinates of disks.
  r
    Disk radius.
  eps
    Minimum spacing between disks and boundary.
  dx
    Bounding box x dimensions.
  dy
    Bounding box y dimensions.

  Returns
  -------

  bool
    Boolean which flags whether all disks edges are above the minimum
    distance from bounding box boundaries.

  """

  # Check bottom, left, top, right bounds
  leftprox = np.abs(c[0]-r) > eps
  rightprox = np.abs(c[0]+r-dx) > eps
  bottomprox = np.abs(c[1]-r) > eps
  topprox = np.abs(c[1]+r-dy) > eps

  # Check corners
  bottomleftprox = np.abs(np.linalg.norm(c)-r) > eps
  br = np.array([dx,0])
  bottomrightprox = np.abs(np.linalg.norm(c-br)-r) > eps
  tr = np.array([dx,dy])
  toprightprox = np.abs(np.linalg.norm(c-tr)-r) > eps
  tl = np.array([0,dy])
  topleftprox = np.abs(np.linalg.norm(c-tl)-r) > eps

  # Return combined check
  return leftprox and rightprox and bottomprox and topprox \
      and bottomleftprox and bottomrightprox and topleftprox and toprightprox



def periodic_copy(m: gmsh.model, c: np.ndarray[float], r: float, dx: float, \
    dy: float, maxtag: int) -> int:
  """Add periodic copies of disks which are over bounding box boundaries.

  Arguments
  ---------

  m
    Gmsh model; shortcut for gmsh.model.
  c
    x and y coordinates of disks.
  r
    Disk radius.
  dx
    Bounding box x dimensions.
  dy
    Bounding box y dimensions.

  Returns
  -------

  int
    Highest Gmsh 2D object tag, after adding periodic disks.


  """

  # Check for disk overlapping boundaries in all 8 periodic copies
  # Add translated disk if so
  if c[0]-r < 0.0:
    m.occ.addDisk(c[0]+dx,c[1],0,r,r,tag=maxtag+1)
    maxtag += 1
  if c[0]+r > dx:
    m.occ.addDisk(c[0]-dx,c[1],0,r,r,tag=maxtag+1)
    maxtag += 1
  if c[1]-r < 0.0:
    m.occ.addDisk(c[0],c[1]+dy,0,r,r,tag=maxtag+1)
    maxtag += 1
  if c[1]+r > dy:
    m.occ.addDisk(c[0],c[1]-dy,0,r,r,tag=maxtag+1)
    maxtag += 1
  if np.linalg.norm(c)-r < 0.0:
    m.occ.addDisk(c[0]+dx,c[1]+dy,0,r,r,tag=maxtag+1)
    maxtag += 1
  if np.linalg.norm(c-np.array([0,dy]))-r < 0.0:
    m.occ.addDisk(c[0]+dx,c[1]-dy,0,r,r,tag=maxtag+1)
    maxtag += 1
  if np.linalg.norm(c-np.array([dx,0]))-r < 0.0:
    m.occ.addDisk(c[0]-dx,c[1]+dy,0,r,r,tag=maxtag+1)
    maxtag += 1
  if np.linalg.norm(c-np.array([dx,dy]))-r < 0.0:
    m.occ.addDisk(c[0]-dx,c[1]-dy,0,r,r,tag=maxtag+1)
    maxtag += 1

  return maxtag



def periodic_disks(nc: int, c: np.ndarray[float], m: gmsh.model, dx: float, \
    dy: float, r: float, eps: float) -> Tuple[Tuple[int,int],int]:
  """Periodic geometry for disks of specified centers and radius.

  Takes circle centers and creates periodically wrapped geometry, with 
  respect to bounding box.

  Arguments
  ---------

  nc
    Number of circles
  m
    Gmsh model; shortcut for gmsh.model.
  c
    x and y coordinates of disks.
  r
    Disk radius.
  dx
    Bounding box x dimensions.
  dy
    Bounding box y dimensions.
  eps
    Minimum spacing between circles.

  Returns
  -------

  boxdimtag
    Tuple containing the dimension and tag number of the bounding box.
  boxtag
    Tag number of the bounding box.

  """

  # Add disks and periodic copies where necessary
  maxtag = 1+nc
  for i,j in enumerate(c):
    m.occ.addDisk(j[0],j[1],0,r,r,tag=i+2)
    maxtag = periodic_copy(m,j,r,dx,dy,maxtag)

  # Fragment overlapping shapes
  out, pc = m.occ.fragment([(2, 1)], [(2, i) for i in range(2, maxtag+1)])

  # Recover the bulk box tag by looking for unique child
  for i in pc[0]:
    bulkbox = True
    for j in pc[1:]:
      for k in j:
        if i == k:
          bulkbox = False
    if bulkbox:
      boxdimtag = i
      break
  boxtag = boxdimtag[1]
  m.occ.synchronize()

  # Eliminate outsiders
  vin = m.getEntitiesInBoundingBox(-eps/2,-eps/2,-eps/2, \
      dx+eps/2,dy+eps/2,eps/2,2)
  for v in vin:
    out.remove(v)
  m.occ.remove(out,True)
  m.occ.synchronize()

  return boxdimtag,boxtag



def periodic_mesh(m,dx,eps):
  """Enforce periodic mesh on x-bounds

  Use Gmsh functionality to make mesh periodically consistent across bounding box
  x-dimension.

  Arguments
  ---------

  m
    Gmsh model; shortcut for gmsh.model.
  dx
    Bounding box x dimensions.
  eps
    Mesh minimum spacing.

  """

  # Loop through boundary lines
  ents = m.getEntities(2)
  boundary = m.getBoundary(ents)
  for i,j in enumerate(boundary):

    # Get start and end coordinates of line
    linepoints = np.array(m.getBoundingBox(j[0],j[1]))
    xmin = linepoints[0]
    xmax = linepoints[3]
    ymin = linepoints[1]

    # Check for vertical line
    if np.abs(xmin-xmax) < eps/2:

      # Loop over all other boundary lines
      for k,l in enumerate(boundary):
        if k != i:
          linepoints2 = np.array(m.getBoundingBox(l[0],l[1]))
          xmin2 = linepoints2[0]
          xmax2 = linepoints2[3]
          ymin2 = linepoints2[1]

          # Check for line directly opposite to set PBCs
          if np.abs(xmin2-xmax2) < eps/2 and \
            np.abs((xmin2-xmin)-dx) < eps/2 and \
            np.abs(ymin-ymin2) < eps/2:
            m.mesh.setPeriodic(1, [j[1]], [l[1]], \
                [1, 0, 0, -dx, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])



def nodeset(f: fi.input, nodes: np.ndarray[int]):
  """Function to write node sets in Abaqus input file

  Arguments
  ---------

  f
    Input file object.
  nodes
    Node set

  """

  # Write node sets to file with max 10 nodes on a line
  ticker = 1
  for i,j in enumerate(nodes):
    if ticker == 10 or i+1 == len(nodes):
      f.write(f'{j},\n')
      ticker = 1
    else:
      f.write(f'{j}, ')
      ticker += 1

