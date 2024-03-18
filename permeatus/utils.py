import numpy as np
import gmsh
import sys
import fileinput as fi
from pandas import unique
import copy
from fractions import Fraction

### Gmsh Functions ###
# Get boundary node sets in 2D box setups
def boundary_nodes_2d(m,dx,dy):

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

# Check that disk outer surface is not too close to a box boundary
# which would distort mesh
def bound_proximity_check_2d(c,r,eps,boxsize):

    # Check bottom, left, top, right bounds
    leftprox = np.abs(c[0]-r) > eps
    rightprox = np.abs(c[0]+r-boxsize) > eps
    bottomprox = np.abs(c[1]-r) > eps
    topprox = np.abs(c[1]+r-boxsize) > eps

    # Check corners
    bottomleftprox = np.abs(np.linalg.norm(c)-r) > eps
    br = np.array([boxsize,0])
    bottomrightprox = np.abs(np.linalg.norm(c-br)-r) > eps
    tr = np.array([boxsize,boxsize])
    toprightprox = np.abs(np.linalg.norm(c-tr)-r) > eps
    tl = np.array([0,boxsize])
    topleftprox = np.abs(np.linalg.norm(c-tl)-r) > eps

    # Return combined check
    return leftprox and rightprox and bottomprox and topprox \
        and bottomleftprox and bottomrightprox and topleftprox and toprightprox

# Add periodic copies of disks
def periodic_copy(m,c,r,boxsize,maxtag):

    # Check for disk overlapping boundaries in all 8 periodic copies
    # Add translated disk if so
    if c[0]-r < 0.0:
        m.addDisk(c[0]+boxsize,c[1],0,r,r,tag=maxtag+1)
        maxtag += 1
    if c[0]+r > boxsize:
        m.addDisk(c[0]-boxsize,c[1],0,r,r,tag=maxtag+1)
        maxtag += 1
    if c[1]-r < 0.0:
        m.addDisk(c[0],c[1]+boxsize,0,r,r,tag=maxtag+1)
        maxtag += 1
    if c[1]+r > boxsize:
        m.addDisk(c[0],c[1]-boxsize,0,r,r,tag=maxtag+1)
        maxtag += 1
    if np.linalg.norm(c)-r < 0.0:
        m.addDisk(c[0]+boxsize,c[1]+boxsize,0,r,r,tag=maxtag+1)
        maxtag += 1
    if np.linalg.norm(c-np.array([0,boxsize]))-r < 0.0:
        m.addDisk(c[0]+boxsize,c[1]-boxsize,0,r,r,tag=maxtag+1)
        maxtag += 1
    if np.linalg.norm(c-np.array([boxsize,0]))-r < 0.0:
        m.addDisk(c[0]-boxsize,c[1]+boxsize,0,r,r,tag=maxtag+1)
        maxtag += 1
    if np.linalg.norm(c-np.array([boxsize,boxsize]))-r < 0.0:
        m.addDisk(c[0]-boxsize,c[1]-boxsize,0,r,r,tag=maxtag+1)
        maxtag += 1
    return maxtag

# Periodic geometry for disks of specified centers and radius
def periodic_disks(nc,centers,m,boxsize,r,eps):

    # Add disks and periodic copies where necessary
    maxtag = 1+nc
    for i,j in enumerate(centers):
        m.occ.addDisk(j[0],j[1],0,r,r,tag=i+2)
        maxtag = periodic_copy(m.occ,j,r,boxsize,maxtag)

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
        boxsize+eps/2,boxsize+eps/2,eps/2,2)
    for v in vin:
        out.remove(v)
    m.occ.remove(out,True)
    m.occ.synchronize()

    return boxdimtag,boxtag

# Enforce periodic mesh on x-bounds
def periodic_mesh(m,boxsize,eps):

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
                        np.abs((xmin2-xmin)-boxsize) < eps/2 and \
                        np.abs(ymin-ymin2) < eps/2:
                        m.mesh.setPeriodic(1, [j[1]], [l[1]], \
                            [1, 0, 0, -boxsize, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])

# Function to write node sets
def nodeset(f,nodes):
    ticker = 1
    for i,j in enumerate(nodes):
        if ticker == 10 or i+1 == len(nodes):
            f.write(f'{j},\n')
            ticker = 1
        else:
            f.write(f'{j}, ')
            ticker += 1

# Abaqus diffusion sim input file production
#TODO add variable inputs for material properties, timepoints, 
# and BC magnitudes
def write_abaqus_diffusion(D,S,C0,C1,touts,tstep,\
    bottomnodes,topnodes,leftnodes,rightnodes,jobname, \
    PBC=True):

    # Replace element types with diffusion
    fname = f'{jobname}.inp'
    gmsh.write(fname)
    with fi.input(fname,inplace=True) as f:
        for line in f:
            print(line.replace("CPS","DC2D"), end='')

    # Append BC's, Diffusion step details, Material properties and section assignment
    with open(fname,"a") as f:
        # Define boundary node sets
        f.write('*NSET, NSET=source \n')
        nodeset(f,bottomnodes)
        f.write('*NSET, NSET=sink \n')
        nodeset(f,topnodes)
        if PBC:
          f.write('*NSET, NSET=lwall, UNSORTED \n')
          nodeset(f,leftnodes)
          f.write('*NSET, NSET=rwall, UNSORTED \n')
          nodeset(f,rightnodes)

        # Periodic concentration for left and right walls
        if PBC:
          f.write('*Equation \n')
          f.write('2 \n')
          f.write("lwall, 11, 1. \n")
          f.write("rwall, 11, -1. \n")

        # Define materials
        for i in range(len(D)):
          f.write(f'*Material, name=material{i}\n')
          f.write(f'*Diffusivity, law=FICK\n')
          f.write(f'{D[i]}, 0.\n')
          f.write(f'*Solubility\n')
          f.write(f'{S[i]},\n')
          f.write(f'*Solid Section, elset=material{i}, material=material{i}\n')

        # Time points
        f.write(f'*Time Points, name=timepoints\n')
        nodeset(f,touts)

        # Zero temperature
        f.write(f'*Physical Constants, absolute zero=0.\n')

        # Diffusion step details
        maxinc = round(touts[-1]/tstep*10)
        f.write(f'*Step, name=diffusion, nlgeom=NO, inc={maxinc}\n')
        f.write(f'*Mass Diffusion, end=PERIOD, dcmax={C0-C1}\n')
        f.write(f'{tstep}, {touts[-1]}, {tstep/10}, {tstep*10},\n')
        f.write(f'*Boundary\n')
        f.write(f'sink, 11, 11, {C1}\n')
        f.write(f'*Boundary\n')
        f.write(f'source, 11, 11, {C0}\n')
        f.write(f'*Restart, write, frequency=0\n')
        f.write(f'*Output, field, time points=timepoints\n')
        f.write(f'*Node Output\n')
        f.write(f'NNC,\n')
        f.write(f'*Element Output, directions=YES\n')
        f.write(f'CONC, MFL, IVOL\n')
        f.write(f'*End Step\n')
