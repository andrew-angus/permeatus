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

# %%
"""
## Gmsh functions
"""

# %%
# Corner fragment
def corner_fragment(m,maxtag,shift,boxsize):
    print(shift,maxtag)
    m.addRectangle(shift[0],shift[1],0,boxsize,boxsize,tag=maxtag+1)
    ents = m.getEntities(2)
    bdimtag = (2,maxtag+1)
    box = [bdimtag]
    ents.remove(bdimtag)
    frags, pc = m.fragment(box, ents)
    bdimtag = pc[0][0]
    box = [bdimtag]
    m.remove(box,recursive=True)
    ents = m.getEntities(2)
    maxtag = np.max([i[1] for i in frags])
    return maxtag

# Periodic shift of overlapping entities
def periodic_shift(m,shift,boxdimtag,eps):
    outs = m.getEntitiesInBoundingBox(-eps/2+shift[0],-eps/2+shift[1],-eps/2,\
                                      boxsize+eps/2+shift[0],boxsize+eps/2+shift[1],eps/2,2)
    m.translate(outs,-shift[0],-shift[1],0)
    ents = m.getEntities(2)
    box = [boxdimtag]
    ents.remove(boxdimtag)
    frags, pc = m.fragment(box, ents)
    boxdimtag = pc[0][0]
    boxtag = boxdimtag[1]
    return frags, boxdimtag, boxtag
    
# Get boundary node sets in 2D box setups
def boundary_nodes_2d(m):
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

#gmsh.option.setNumber("Mesh.SaveGroupsOfElements", 1)
gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)
#gmsh.option.setNumber("Geometry.OCCBooleanPreserveNumbering", 1)

# Bounding box
boxsize = 1
gmsh.model.occ.addRectangle(0,0,0,boxsize,boxsize,tag=1)

# Circel radius and inimum spacing
r = 0.2
eps = r/10

# Periodic translation arrays
translationbase = np.array([-boxsize,0.0,boxsize])
translations = np.array([[i,j] for i in translationbase for j in translationbase])
translations = np.delete(translations,4,axis=0)
corners = np.array([[i,j] for i in translationbase[::2] for j in translationbase[::2]])

# Add a randomly centred circle
while True:
    #c = np.random.rand(2)*np.array([boxsize,boxsize]) # Random position
    c = np.array([0.0,0.0]) # Corner test case
    # Check for minimum distance between circle exterior and bounding box
    if np.abs(c[0]-r) > eps and np.abs(c[0]+r-boxsize) > eps and \
        np.abs(c[1]-r) > eps and np.abs(c[1]+r-boxsize) > eps:
        gmsh.model.occ.addDisk(c[0],c[1],0,r,r,tag=2)
        break

# Fragment overlapping shapes
frags, pc = gmsh.model.occ.fragment([(2, 1)], [(2, i) for i in range(2, 3)])
boxdimtag = pc[0][0]
boxtag = boxdimtag[1]
maxtag = np.max([i[1] for i in frags])

# Also frament corner periodic translations
for i in corners:
    maxtag = corner_fragment(gmsh.model.occ,maxtag,i,boxsize)

# Periodic shifts
for i in translations:
    frags,boxdimtag, boxtag = periodic_shift(gmsh.model.occ,i,boxdimtag,eps)

# Synchronise
gmsh.model.occ.synchronize()

# Identify physical groups for material assignment
frags.remove(boxdimtag)
gmsh.model.addPhysicalGroup(2, [boxtag], name="layer0")
gmsh.model.addPhysicalGroup(2, [i[1] for i in frags], name="layer1")

# Generate mesh
gmsh.model.mesh.generate(2)

# Acquire boundary node sets
gmsh.model.occ.synchronize()
bottomnodes, topnodes, leftnodes, rightnodes = boundary_nodes_2d(gmsh.model)

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

#gmsh.option.setNumber("Mesh.SaveGroupsOfElements", 1)
gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)
#gmsh.option.setNumber("Geometry.OCCBooleanPreserveNumbering", 1)

# Microstructure specifications
vfrac = 0.5
r = 5.2e-6
area0 = np.pi*r**2
nc = 15
area1 = nc*area0
area2 = area1/vfrac

# Bounding box
boxsize = np.sqrt(area2)
print(boxsize*1e6)
gmsh.model.occ.addRectangle(0,0,0,boxsize,boxsize,tag=1)

# Minimum spacing
eps = r*0.05
cbuff = 2*r + eps

# Translation arrays
translationbase = np.array([-boxsize,0.0,boxsize])
translations = np.array([[i,j] for i in translationbase for j in translationbase])
translations = np.delete(translations,4,axis=0)
corners = np.array([[i,j] for i in translationbase[::2] for j in translationbase[::2]])

# Add randomly centred circles
# Loop over number of desired circles
centers = np.empty((0,2))
for i in range(nc):

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
                    mindist = np.linalg.norm(c[0]-center)
                    for translator in translations:
                        dist = np.linalg.norm(c[0]+translator-center)
                        mindist = np.minimum(dist,mindist)
                    if mindist < cbuff:
                        reject = True
        niter += 1
        if niter > 100000:
            raise Exception("Structure generation took too many iterations")

    # Add accepted circle
    print(i,c)
    centers = np.r_[centers,c]
    gmsh.model.occ.addDisk(c[0,0],c[0,1],0,r,r,tag=2+i)

# Fragment overlapping shapes
frags, pc = gmsh.model.occ.fragment([(2, 1)], [(2, i) for i in range(2, 2+nc)])
boxdimtag = pc[0][0]
boxtag = boxdimtag[1]
maxtag = np.max([i[1] for i in frags])

# Also frament corner periodic translations
for i in corners:
    maxtag = corner_fragment(gmsh.model.occ,maxtag,i,boxsize)

# Periodic shifts
for i in translations:
    frags,boxdimtag, boxtag = periodic_shift(gmsh.model.occ,i,boxdimtag,eps)

# Synchronise
gmsh.model.occ.synchronize()

#TODO Need to fix tag references to account for fragmenting
frags.remove(boxdimtag)
gmsh.model.addPhysicalGroup(2, [boxtag], name="layer0")
gmsh.model.addPhysicalGroup(2, [i[1] for i in frags], name="layer1")

# We can then generate a 2D mesh...
gmsh.option.setNumber("Mesh.MeshSizeMin",eps)
gmsh.option.setNumber("Mesh.MeshSizeMax",r/5)
gmsh.model.mesh.generate(2)

gmsh.model.occ.synchronize()
bottomnodes, topnodes, leftnodes, rightnodes = boundary_nodes_2d(gmsh.model)

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

#gmsh.option.setNumber("Mesh.SaveGroupsOfElements", 1)
gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)
#gmsh.option.setNumber("Geometry.OCCBooleanPreserveNumbering", 1)

# Microstructure specifications
#minspacefrac = 1e-3
vfrac = 0.5
r = 0.03
area0 = np.pi*r**2
nc = 100
area1 = nc*area0
area2 = area1/vfrac

# Bounding box
boxsize = np.sqrt(area2)
print(boxsize)
gmsh.model.occ.addRectangle(0,0,0,boxsize,boxsize,tag=1)

# Minimum spacing
eps = 1e-2*boxsize

# Define buffers
ybuff = r + eps
cbuff = 2*r + eps

# Run LS algorithm
# Randomly insert N point particles
c = np.random.rand(nc,2)*np.array([[boxsize,boxsize-2*ybuff]])
c[:,1] += ybuff

# Randomly assign velocities
v = np.random.rand(nc,2)*np.array([[2*boxsize,2*(boxsize-2*ybuff)]])
v[:,0] -= boxsize
v[:,1] -= boxsize-2*ybuff

plt.scatter(c[:,0],c[:,1])
plt.show()

# Define algorithm variables
h = 1e-2 # Radius growth rate
t = 0 # Time
rt = 0 # r(t)
ymax = boxsize - eps/2 # Max y-coordinate
ymin = eps/2 # Min y-coordinate
rf = r + eps/2 # Final radius, including minimum spacing
tf = rf/h # Final time
wallhits = np.ones(nc)*np.inf
parthits = np.ones((nc,nc))*np.inf
cp = np.zeros((nc,2))
print(ymin,ymax)
niter = 0

# Periodic wrapping of x-coordinate
def xwrap(point):
    point[0] -= np.floor(point[0]/boxsize)*boxsize
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
while advancing and niter < 1000000:

    # List of wall collisions
    for i in range(nc):
        if v[i,1] < 0.0:
            bound = ymin
        else:
            bound = ymax
        wallhits[i] = (bound - np.sign(v[i,1])*rt - c[i,1])/(v[i,1]+np.sign(v[i,1])*h)
        cp[i] = xwrap(c[i] + v[i]*wallhits[i] + np.sign(v[i,1])*np.array([0,rt+h*wallhits[i]]))

    # Find minimum
    bouncer = np.argmin(wallhits)
    time = wallhits[bouncer]
    wallhit = True

    # Get particle-particle collisions
    for i in range(nc):
        for j in range(i):
            # Get central time
            tcoll = ppColl(c[i],c[j],v[i],v[j],h,rt)
            tcoll = np.minimum(ppColl(c[i]-np.array([boxsize,0]),c[j],v[i],v[j],h,rt),tcoll)
            parthits[i,j] = np.minimum(ppColl(c[i],c[j]-np.array([boxsize,0]),v[i],v[j],h,rt),tcoll)
            if parthits[i,j] < 0:
                print(f'Times: {parthits[i,j]}, {ppColl(c[i],c[j],v[i],v[j],h,rt)}, \
                {ppColl(c[i]-np.array([boxsize,0]),c[j],v[i],v[j],h,rt)}, {ppColl(c[i],c[j]-np.array([boxsize,0]),v[i],v[j],h,rt)}')
                print(f'Particles: {i}, {j}')
                print(f'Original positions: {c[i]}, {c[j]}')
                # Propagate centers
                c[i] = xwrap(c[i] + v[i]*time)
                c[j] = xwrap(c[j] + v[j]*time)
                dc = c[colliders[0]] - c[colliders[1]]
                dc[0] -= round(dc[0] / boxsize) * boxsize
                u = dc/np.linalg.norm(dc)
                print(f'Velocities: {v[i]}, {v[j]}')
                print(f'Final positions: {c[i]}, {c[j]}')
                print(f'Impact point: {(c[i]+c[j])/2}')
                print(f'Impact normal: {u}')
                vpi = np.dot(u,v[i])*u
                vti = v[i] - vpi
                vpj = np.dot(u,v[j])*u
                vtj = v[j] - vpj
                print(f'Parallel velocities: {vpi}, {vpj}')
                print(f'Tangential velocities: {vti}, {vtj}')
                v[i] = vpj + vti
                v[j] = vpi + vtj
                print(f'Pre-boosts: {v[i]}, {v[j]}')
                v[i] += 2*u*h
                v[j] -= 2*u*h
                print(f'Finals: {v[i]}, {v[j]}')
                
                print('')
                break

    # Find minimum and compare with wallhit time
    colliders = np.unravel_index(parthits.argmin(), parthits.shape)
    pptime = parthits[colliders]
    if pptime < time:
        wallhit = False
        time = pptime


    # Check for final time
    if t + time > tf:
        time = tf - t 
        for i in range(nc):
            c[i] = xwrap(c[i]+v[i]*time)
        advancing = False
    else:
        if wallhit:
            print(f'Particle: {bouncer}')
            print(f'Original position: {c[bouncer]}')
            c[bouncer] = xwrap(c[bouncer] + v[bouncer]*time) # Propagate c
            u = (c[bouncer]-cp[bouncer])/np.linalg.norm(c[bouncer]-cp[bouncer])
            print(f'Velocity: {v[bouncer]}')
            print(f'Propagation time: {time}')
            print(f'Final position: {c[bouncer]}')
            print(f'Impact point: {cp[bouncer]}')
            print(f'Impact normal: {u}')
            print(f'Velocity normal: {v[bouncer]/np.linalg.norm(v[bouncer])}')
            vp = np.dot(u,v[bouncer])*u
            vt = v[bouncer] - vp
            print(f'Parallel velocity: {vp}')
            print(f'Tangential velocity: {vt}')
            v[bouncer,1] *= -1
            print(f'Pre-boost: {v[bouncer]}')
            v[bouncer] += u*h
            print(f'Final velocity: {v[bouncer]}')
            print('')
        
        else:
            print(f'Particles: {colliders[0]}, {colliders[1]}')
            print(f'Propagation time: {time}')
            print(f'Start and final radius: {rt}, {rt + time*h}')
            print(f'Original positions: {c[colliders[0]]}, {c[colliders[1]]}')
            # Propagate centers
            c[colliders[0]] = xwrap(c[colliders[0]] + v[colliders[0]]*time)
            c[colliders[1]] = xwrap(c[colliders[1]] + v[colliders[1]]*time)
            dc = c[colliders[0]] - c[colliders[1]]
            dc[0] -= round(dc[0] / boxsize) * boxsize
            u = dc/np.linalg.norm(dc)
            print(f'Velocities: {v[colliders[0]]}, {v[colliders[1]]}')
            print(f'Final positions: {c[colliders[0]]}, {c[colliders[1]]}')
            print(f'Impact point: {(c[colliders[0]]+c[colliders[1]])/2}')
            print(f'Impact normal: {u}')
            vpi = np.dot(u,v[colliders[0]])*u
            vti = v[colliders[0]] - vpi
            vpj = np.dot(u,v[colliders[1]])*u
            vtj = v[colliders[1]] - vpj
            print(f'Parallel velocities: {vpi}, {vpj}')
            print(f'Tangential velocities: {vti}, {vtj}')
            v[colliders[0]] = vpj + vti
            v[colliders[1]] = vpi + vtj
            print(f'Pre-boosts: {v[colliders[0]]}, {v[colliders[1]]}')
            v[colliders[0]] += 2*u*h
            v[colliders[1]] -= 2*u*h
            print(f'Finals: {v[colliders[0]]}, {v[colliders[1]]}')
            print('')
            
    if time < 0:
        print('exiting ', niter, t, rt)
        break
    
    # Propagate non-bouncing centers
    for i in range(nc):
        if (wallhit and i != bouncer) or ((not wallhit) and i != colliders[0] and i != colliders[1]):
            c[i] = xwrap(c[i]+v[i]*time)
            
    # Update trackers
    t += time
    rt += h*time
    niter += 1

    print(niter,t,rt)
    print('')

    # Cap time to ybuff time, else next collision

plt.scatter(c[:,0],c[:,1])
plt.show()

for i in range(nc):
    #gmsh.model.occ.addDisk(c[i,0],c[i,1],0,r,r,tag=2+i)
    gmsh.model.occ.addDisk(c[i,0],c[i,1],0,r,r,tag=2+i)

# Fragment overlapping shapes
frags, pc = gmsh.model.occ.fragment([(2, 1)], [(2, i) for i in range(2, 2+nc)])
boxdimtag = pc[0][0]
boxtag = boxdimtag[1]
maxtag = np.max([i[1] for i in frags])

# Also frament corner periodic translations
for i in corners:
    maxtag = corner_fragment(gmsh.model.occ,maxtag,i,boxsize)

# Periodic shifts
for i in translations:
    frags,boxdimtag, boxtag = periodic_shift(gmsh.model.occ,i,boxdimtag,eps)

# Synchronise
gmsh.model.occ.synchronize()

#TODO Need to fix tag references to account for fragmenting
frags.remove(boxdimtag)
gmsh.model.addPhysicalGroup(2, [boxtag], name="layer0")
#gmsh.model.removeEntities([boxdimtag],True)
gmsh.model.addPhysicalGroup(2, [i[1] for i in frags], name="layer1")

# We can then generate a 2D mesh...
gmsh.option.setNumber("Mesh.MeshSizeMin",eps)
gmsh.option.setNumber("Mesh.MeshSizeMax",r/5)
gmsh.model.mesh.generate(2)

gmsh.model.occ.synchronize()
bottomnodes, topnodes, leftnodes, rightnodes = boundary_nodes_2d(gmsh.model)

gmsh.write("gmsh.msh")
write_abaqus_diffusion(bottomnodes, topnodes, leftnodes, rightnodes)
gmsh.fltk.run()
gmsh.finalize()

# %%
0.001126278067159134