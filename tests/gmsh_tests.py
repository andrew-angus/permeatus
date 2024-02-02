# %%
"""
# Setup
"""

# %%
import numpy as np
import gmsh
import sys
import fileinput as fi
from pandas import unique
import matplotlib.pyplot as plt
from tqdm.notebook import trange
import copy

# %%
"""
## Gmsh functions
"""

# %%
# Get boundary node sets in 2D box setups
def boundary_nodes_2d(m,boxsize):
    boundary = m.getBoundary(gmsh.model.getEntities(2))
    leftnodes = np.empty(0)
    leftcoords = np.empty(0)
    rightnodes = np.empty(0)
    rightcoords = np.empty(0)
    topnodes = np.empty(0)
    topcoords = np.empty(0)
    bottomnodes = np.empty(0)
    bottomcoords = np.empty(0)
    for dim, tag in boundary:
        nodeTags, coord, parametricCoord = m.mesh.getNodes(dim, np.abs(tag), includeBoundary=True)
        coords = np.reshape(coord, (len(coord)//3, 3))
        main_coord = np.argmax(np.abs(np.sum(np.diff(coords,axis=0),axis=0)))
        if main_coord == 0:
            if np.isclose(coord[1],0.0):
                bottomnodes = np.r_[bottomnodes,nodeTags]
                bottomcoords = np.r_[bottomcoords, coords[:,main_coord]]
            elif np.isclose(coord[1],boxsize):
                topnodes = np.r_[topnodes,nodeTags]
                topcoords = np.r_[topcoords, coords[:,main_coord]]
        elif main_coord == 1:
            if np.isclose(coord[0],0.0):
                leftnodes = np.r_[leftnodes,nodeTags]
                leftcoords = np.r_[leftcoords, coords[:,main_coord]]
            elif np.isclose(coord[0],boxsize):
                rightnodes = np.r_[rightnodes,nodeTags]
                rightcoords = np.r_[rightcoords, coords[:,main_coord]]
    
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

    return bottomnodes, topnodes, leftnodes, rightnodes

def bound_proximity_check_2d(c,r,eps,boxsize):
    leftprox = np.abs(c[0]-r) > eps
    rightprox = np.abs(c[0]+r-boxsize) > eps
    bottomprox = np.abs(c[1]-r) > eps
    topprox = np.abs(c[1]+r-boxsize) > eps
    bottomleftprox = np.abs(np.linalg.norm(c)-r) > eps
    br = np.array([boxsize,0])
    bottomrightprox = np.abs(np.linalg.norm(c-br)-r) > eps
    tr = np.array([boxsize,boxsize])
    toprightprox = np.abs(np.linalg.norm(c-tr)-r) > eps
    tl = np.array([0,boxsize])
    topleftprox = np.abs(np.linalg.norm(c-tl)-r) > eps
    return leftprox and rightprox and bottomprox and topprox \
        and bottomleftprox and bottomrightprox and topleftprox and toprightprox
    
def periodic_copy(m,c,r,boxsize,maxtag):
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
    out, pc = gmsh.model.occ.fragment([(2, 1)], [(2, i) for i in range(2, maxtag+1)])
    boxdimtag = pc[0][0]
    boxtag = boxdimtag[1]
    m.occ.synchronize()
    
    # Eliminate outsiders
    vin = m.getEntitiesInBoundingBox(-eps/2,-eps/2,-eps/2,boxsize+eps/2,boxsize+eps/2,eps/2,2)
    for v in vin:
        out.remove(v)
    m.occ.remove(out,True)
    m.occ.synchronize()
    
    return boxdimtag,boxtag

# Enforce periodic mesh on x-bounds
def periodic_mesh(m,boxsize,eps):
    
    ents = m.getEntities(2)
    boundary = m.getBoundary(ents)
    for i,j in enumerate(boundary):
        linepoints = np.array(m.getBoundingBox(j[0],j[1]))
        xmin = linepoints[0]
        xmax = linepoints[3]
        ymin = linepoints[1]
        if np.abs(xmin-xmax) < eps/2:
            for k,l in enumerate(boundary):
                if k != i:
                    linepoints2 = np.array(m.getBoundingBox(l[0],l[1]))
                    xmin2 = linepoints2[0]
                    xmax2 = linepoints2[3]
                    ymin2 = linepoints2[1]
                    if np.abs(xmin2-xmax2) < eps/2 and np.abs((xmin2-xmin)-boxsize) < eps/2 and np.abs(ymin-ymin2) < eps/2:
                        m.mesh.setPeriodic(1, [j[1]], [l[1]], [1, 0, 0, -boxsize, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
                        

# %%
"""
## Abaqus Output File Functions
"""

# %%
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
#TODO add variable inputs for material properties, timepoints, and BC magnitudes (best integrating with permeatus)
def write_abaqus_diffusion(sourcenodes,sinknodes,leftnodes,rightnodes,fname="gmsh.inp"):
    # Replace element types with diffusion
    gmsh.write(fname)
    with fi.input(fname,inplace=True) as f:
        for line in f:
            print(line.replace("CPS","DC2D"), end='')
    
    # Append BC's, Diffusion step details, Material properties and section assignment
    with open(fname,"a") as f:
        # Define boundary node sets
        f.write('*NSET, NSET=source \n')
        nodeset(f,sourcenodes)
        f.write('*NSET, NSET=sink \n')
        nodeset(f,sinknodes)
        f.write('*NSET, NSET=lwall, UNSORTED \n')
        nodeset(f,leftnodes)
        f.write('*NSET, NSET=rwall, UNSORTED \n')
        nodeset(f,rightnodes)

        # Periodic concentration for left and right walls
        f.write('*Equation \n')
        f.write('2 \n')
        f.write("lwall, 11, 1. \n")
        f.write("rwall, 11, -1. \n")

    
        # Define materials
        f.write(f'*Material, name=layer0\n')
        f.write(f'*Diffusivity, law=FICK\n')
        f.write(f'1.,0.\n')
        f.write(f'*Solubility\n')
        f.write(f'1.,\n')
        f.write(f'*Material, name=layer1\n')
        f.write(f'*Diffusivity, law=FICK\n')
        f.write(f'0.1,0.\n')
        f.write(f'*Solubility\n')
        f.write(f'1.1,\n')
    
        # Assign material sections
        f.write(f'*Solid Section, elset=layer0, material=layer0\n')
        f.write(f'*Solid Section, elset=layer1, material=layer1\n')
    
        # Time points
        f.write(f'*Time Points, name=timepoints\n')
        f.write(f'0.001, 0.05, 0.2, 2.\n')
    
        # Zero temperature
        f.write(f'*Physical Constants, absolute zero=0.\n')
    
        # Diffusion step details
        f.write(f'*Step, name=diffusion, nlgeom=NO, inc=200000\n')
        f.write(f'*Mass Diffusion, end=PERIOD, dcmax=1.\n')
        f.write(f'0.001, 2., 0.000175, 0.001,\n')
        f.write(f'*Boundary\n')
        f.write(f'sink, 11, 11\n')
        f.write(f'*Boundary\n')
        f.write(f'source, 11, 11, 1.\n')
        f.write(f'*Restart, write, frequency=0\n')
        f.write(f'*Output, field, time points=timepoints\n')
        f.write(f'*Element Output, position=NODES, directions=YES\n')
        f.write(f'CONC, MFL\n')
        f.write(f'*End Step\n')
    
    # Check
    #with open(fname,"r") as f:
    #    print(f.read())

# %%
"""
# ABAQUS 2-layer
"""

# %%
# Initialise
gmsh.initialize()
gmsh.model.add("gmsh")

# Output both element and node groups
gmsh.option.setNumber("Mesh.SaveGroupsOfElements", 1)
gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)

# Geometry points
width = 1.0/80
gmsh.model.geo.addPoint(0, 0, 0, tag=1,meshSize=width)
gmsh.model.geo.addPoint(width*10, 0, 0, tag=2)#,meshSize=width)
gmsh.model.geo.addPoint(width*10, .5, 0, tag=3)#,meshSize=width)
gmsh.model.geo.addPoint(0, .5, 0, tag=4)#,meshSize=width)
gmsh.model.geo.addPoint(width*10, 1.0, 0, tag=5)#,meshSize=width)
gmsh.model.geo.addPoint(.0, 1.0, 0, tag=6)#,meshSize=width)

# Lines
gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(2, 3, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)
gmsh.model.geo.addLine(3, 5, 5)
gmsh.model.geo.addLine(5, 6, 6)
gmsh.model.geo.addLine(6, 4, 7)

# Curve loops - Need to both be anti-clockwise, or both clockwise
gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
gmsh.model.geo.addCurveLoop([-3, 5, 6, 7], 2)

# Add surfacets
gmsh.model.geo.addPlaneSurface([1], 1)
gmsh.model.geo.addPlaneSurface([2], 2)

# Define structured mesh with seeded edges
gmsh.model.geo.mesh.setTransfiniteCurve(1, 2, "Progression", 1.0)
gmsh.model.geo.mesh.setTransfiniteCurve(6, 2)
gmsh.model.geo.mesh.setTransfiniteCurve(3, 2)
gmsh.model.geo.mesh.setTransfiniteCurve(2, 40)
gmsh.model.geo.mesh.setTransfiniteCurve(4, 40)
gmsh.model.geo.mesh.setTransfiniteCurve(5, 40)
gmsh.model.geo.mesh.setTransfiniteCurve(7, 40)
gmsh.model.geo.mesh.setTransfiniteSurface(1)
gmsh.model.geo.mesh.setTransfiniteSurface(2)

# Finalise geometry
gmsh.model.geo.synchronize()

# Define physical groups
gmsh.model.addPhysicalGroup(2, [1], name="layer0")
gmsh.model.addPhysicalGroup(2, [2], name="layer1")

# Generate mesh
#gmsh.option.setNumber("Mesh.MeshSizeMin",width)
#gmsh.option.setNumber("Mesh.MeshSizeMax",width)
#gmsh.model.mesh.setAlgorithm(2,2,8)
gmsh.model.mesh.generate(2)
gmsh.model.mesh.recombine()

# Obtain source and sink nodes
sourcenodes = gmsh.model.mesh.getElements(1,1)[-1][0]
sinknodes = gmsh.model.mesh.getElements(1,6)[-1][0]
leftnodes = np.append(gmsh.model.mesh.getElements(1,7)[-1][0], gmsh.model.mesh.getElements(1,4)[-1][0]) 
rightnodes = np.append(gmsh.model.mesh.getElements(1,2)[-1][0], gmsh.model.mesh.getElements(1,5)[-1][0])
print(sourcenodes)
print(sinknodes)
print(leftnodes)
leftnodes = np.array([i for i in leftnodes if i not in sourcenodes])
leftnodes = np.array([i for i in leftnodes if i not in sinknodes])
print(leftnodes)
print(rightnodes)
rightnodes = np.array([i for i in rightnodes if i not in sourcenodes])
rightnodes = np.array([i for i in rightnodes if i not in sinknodes])
print(rightnodes)
sourcenodes = unique(sourcenodes)
sinknodes = unique(sinknodes)
leftnodes = unique(leftnodes)
rightnodes = unique(np.flip(rightnodes))
print(sourcenodes)
print(sinknodes)
print(leftnodes)
print(rightnodes)
#print(gmsh.model.mesh.getElements())

# Save for opening in gmsh
gmsh.write("gmsh2layer.msh")

# Save for running ABAQUS sim
write_abaqus_diffusion(sourcenodes,sinknodes,leftnodes,rightnodes)

# Visualise
gmsh.fltk.run()

# Finalise
gmsh.finalize()

# %%
"""
# Box with one central hole
"""

# %%
"""
# Box with random hole
"""

# %%
try:
    gmsh.finalize()
except:
    pass
gmsh.initialize()

gmsh.model.add("random")

gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)

# Bounding box
boxsize = 1
gmsh.model.occ.addRectangle(0,0,0,boxsize,boxsize,tag=1)

# Circel radius and inimum spacing
r = 0.2
eps = r/20

# Translations
translationbase = np.array([-boxsize,0.0,boxsize])
translations = np.array([[i,j] for i in translationbase for j in translationbase])

# Test case override
#c = np.array([[0.0,0.0],[2*r+eps,0],[0,boxsize/2]])
c = np.array([[0,0]])
nc = len(c)

# Add a randomly centred circle
"""
while True:
    #c = np.random.rand(2)*np.array([boxsize,boxsize]) # Random position
    c = np.array([0.0,0.0]) # Corner test case
    # Check for minimum distance between circle exterior and bounding box
    if np.abs(c[0]-r) > eps and np.abs(c[0]+r-boxsize) > eps and \
        np.abs(c[1]-r) > eps and np.abs(c[1]+r-boxsize) > eps:
        gmsh.model.occ.addDisk(c[0],c[1],0,r,r,tag=2)
        break
"""

# Add circles with periodic wrapping
boxdimtag,boxtag = periodic_disks(nc,c,gmsh.model,boxsize,r,eps)

# Identify physical groups for material assignment
ents = gmsh.model.getEntities(2)
ents.remove(boxdimtag)
gmsh.model.addPhysicalGroup(2, [boxtag], name="layer0")
gmsh.model.addPhysicalGroup(2, [i[1] for i in ents], name="layer1")

# Enforce periodic mesh on x-bounds
periodic_mesh(gmsh.model,boxsize,eps)

# Generate mesh
gmsh.model.mesh.generate(2)

# Acquire boundary node sets
gmsh.model.occ.synchronize()
bottomnodes, topnodes, leftnodes, rightnodes = boundary_nodes_2d(gmsh.model,boxsize)
print(leftnodes,rightnodes)

gmsh.write("gmsh.msh")
write_abaqus_diffusion(bottomnodes,topnodes,leftnodes,rightnodes)
gmsh.fltk.run()
gmsh.finalize()

# %%
"""
# Box with more random holes
"""

# %%
try:
    gmsh.finalize()
except:
    pass
gmsh.initialize()

gmsh.model.add("random")

gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)
#gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)

# Microstructure specifications
vfrac = 0.4
r = 0.1
area0 = np.pi*r**2
nc = 15
area1 = nc*area0
area2 = area1/vfrac

# Bounding box
boxsize = np.sqrt(area2)
print(boxsize)
gmsh.model.occ.addRectangle(0,0,0,boxsize,boxsize,tag=1)

# Minimum spacing
eps = r/20
cbuff = 2*r + eps

# Translations
translationbase = np.array([-boxsize,0.0,boxsize])
translations = np.array([[i,j] for i in translationbase for j in translationbase])

# Add randomly centred circles
# Loop over number of desired circles
centers = np.empty((0,2))
for i in trange(nc):

    # Continually try and insert circle within contstraints
    niter = 0
    reject = True
    while reject:

        # Randomly draw circle center
        c = np.random.rand(1,2)*np.array([[boxsize,boxsize]])
        
        # For first circle just check minimum distance from edges and corners
        if bound_proximity_check_2d(c[0],r,eps,boxsize):
            if i == 0:
                reject = False
                        
            # Afterwards check against distance between previous circles
            else:
                reject = False
                for center in centers:
                    mindist = np.inf
                    for translator in translations:
                        dist = np.linalg.norm(c[0]+translator-center)
                        mindist = np.minimum(dist,mindist)
                    if mindist < cbuff:
                        reject = True
        niter += 1
        if niter > 5000:
            raise Exception("Structure generation took too many iterations")

    # Add accepted circle
    centers = np.r_[centers,c]

print(niter)

# Add circles with periodic wrapping
boxdimtag,boxtag = periodic_disks(nc,centers,gmsh.model,boxsize,r,eps)

# Identify physical groups for material assignment
ents = gmsh.model.getEntities(2)
ents.remove(boxdimtag)
gmsh.model.addPhysicalGroup(2, [boxtag], name="layer0")
gmsh.model.addPhysicalGroup(2, [i[1] for i in ents], name="layer1")

# Enforce periodic mesh on x-bounds
periodic_mesh(gmsh.model,boxsize,eps)

# Scale mesh to microns
scaling = 5.2e-6/0.1
gmsh.model.occ.dilate(gmsh.model.getEntities(2),0,0,0,scaling,scaling,scaling)
gmsh.model.occ.synchronize()

# Generate mesh
gmsh.option.setNumber("Mesh.MeshSizeMin",eps*scaling)
gmsh.option.setNumber("Mesh.MeshSizeMax",r*scaling/5)
gmsh.model.mesh.generate(2)

# Acquire boundary node sets
gmsh.model.occ.synchronize()
bottomnodes, topnodes, leftnodes, rightnodes = boundary_nodes_2d(gmsh.model,boxsize*scaling)

gmsh.write("gmsh.msh")
write_abaqus_diffusion(bottomnodes,topnodes,leftnodes,rightnodes)
gmsh.fltk.run()
gmsh.finalize()

# %%
"""
# LS Algorithm
"""

# %%
try:
    gmsh.finalize()
except:
    pass
gmsh.initialize()

gmsh.model.add("microstructure")

gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)
#gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)

# Microstructure specifications
vfrac = 0.8
r = 0.1
area0 = np.pi*r**2
nc = 15
area1 = nc*area0
area2 = area1/vfrac

# Bounding box
boxsize = np.sqrt(area2)
print(boxsize)
gmsh.model.occ.addRectangle(0,0,0,boxsize,boxsize,tag=1)

# Minimum spacing
eps = r/20

# Define buffers
cbuff = 2*r + eps

# Run LS algorithm
# Randomly insert N point particles
c = np.random.rand(nc,2)*np.array([[boxsize,boxsize]])

# Randomly assign velocities
v = np.random.rand(nc,2)*np.array([[2*boxsize,2*boxsize]])
v[:,0] -= boxsize
v[:,1] -= boxsize

# Define algorithm variables
h = r # Radius growth rate
t = 0 # Time
rt = 0 # r(t)
rf = r + eps/2 # Final radius, including minimum spacing
tf = rf/h # Final time
parthits = np.ones((nc,nc))*np.inf
niter = 0

# Output times
tout = np.linspace(1e-1,tf,10)
radii = tout*h
output = 0
cout = {}

# Translations
translationbase = np.array([-boxsize,0.0,boxsize])
translations = np.array([[i,j] for i in translationbase for j in translationbase])

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

    # Get particle-particle collisions
    for i in range(nc):
        for j in range(i):
            # Get collision times with each periodic translation and take minimum
            parthits[i,j] = np.inf
            for translator in translations:
                tcoll = ppColl(c[i]+translator,c[j],v[i],v[j],h,rt)
                parthits[i,j] = np.minimum(parthits[i,j],tcoll)

    # Find minimum
    colliders = np.unravel_index(parthits.argmin(), parthits.shape)
    time = parthits[colliders]
    
    # Check for output time
    if t + time > tout[output]:
        time = tout[output] - t 
        for i in range(nc):
            c[i] = wrap(c[i]+v[i]*time)
        cout[output] = copy.deepcopy(c)
        output += 1
        print(output,niter,t+time,rt+h*time)
        if output == len(tout):
            advancing = False
        
        # Plot status
        fig, ax = plt.subplots()
        for i in c:
            circle = plt.Circle((i[0], i[1]), rt+h*time,label=f't = {t+time} s')
            ax.add_patch(circle)
        ax.set_box_aspect(1)
        plt.xlim(0,boxsize)
        plt.ylim(0,boxsize)
        plt.show()
    else:
        # Propagate centers
        c[colliders[0]] = wrap(c[colliders[0]] + v[colliders[0]]*time)
        c[colliders[1]] = wrap(c[colliders[1]] + v[colliders[1]]*time)
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
            
        if time < 0:
            print('exiting ', niter, t, rt)
            break
        
        # Propagate non-colliding centers
        for i in range(nc):
            if i != colliders[0] and i != colliders[1]:
                c[i] = wrap(c[i]+v[i]*time)
            
    # Update trackers
    t += time
    rt += h*time
    niter += 1

    # Check
    if niter > 10000:
        print(niter,t,rt)
        raise Exception("Structure generation took too many iterations")

print(niter,t,rt)

# Add circles with periodic wrapping
boxdimtag,boxtag = periodic_disks(nc,c,gmsh.model,boxsize,r,eps)

# Identify physical groups for material assignment
ents = gmsh.model.getEntities(2)
ents.remove(boxdimtag)
gmsh.model.addPhysicalGroup(2, [boxtag], name="layer0")
gmsh.model.addPhysicalGroup(2, [i[1] for i in ents], name="layer1")

# Enforce periodic mesh on x-bounds
periodic_mesh(gmsh.model,boxsize,eps)

# Scale mesh to microns
scaling = 5.2e-6/r
gmsh.model.occ.dilate(gmsh.model.getEntities(2),0,0,0,scaling,scaling,scaling)
gmsh.model.occ.synchronize()

# Generate mesh
gmsh.option.setNumber("Mesh.MeshSizeMin",eps*scaling)
gmsh.option.setNumber("Mesh.MeshSizeMax",r*scaling/5)
gmsh.model.mesh.generate(2)

# Acquire boundary node sets
gmsh.model.occ.synchronize()
bottomnodes, topnodes, leftnodes, rightnodes = boundary_nodes_2d(gmsh.model,boxsize*scaling)

gmsh.write("gmsh.msh")
write_abaqus_diffusion(bottomnodes,topnodes,leftnodes,rightnodes)
gmsh.fltk.run()
gmsh.finalize()

# %%
