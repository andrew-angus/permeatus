# %%
"""
# ABAQUS 2-layer
"""

# %%
import gmsh
import sys
import fileinput as fi
import numpy as np

# Abaqus diffusion sim input file production
#TODO add variable inputs for material properties, timepoints, and BC magnitudes (best integrating with permeatus)
def write_abaqus_diffusion(sourcenodes,sinknodes,fname="gmsh.inp"):
    # Replace element types with diffusion
    gmsh.write(fname)
    with fi.input(fname,inplace=True) as f:
        for line in f:
            print(line.replace("CPS","DC2D"), end='')
    
    # Append BC's, Diffusion step details, Material properties and section assignment
    with open(fname,"a") as f:
        # Define boundary node sets
        f.write('*NSET, NSET=source \n')
        ticker = 1
        for j in sourcenodes:
            if ticker == 10 or ticker == len(sourcenodes):
                f.write(f'{j},\n')
                ticker = 1
            else:
                f.write(f'{j}, ')
                ticker += 1
        f.write('*NSET, NSET=sink \n')
        ticker = 1
        for j in sinknodes:
            if ticker == 10 or ticker == len(sinknodes):
                f.write(f'{j},\n')
                ticker = 1
            else:
                f.write(f'{j}, ')
                ticker += 1
    
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
#print(gmsh.model.mesh.getElements())

# Save for opening in gmsh
gmsh.write("gmsh2layer.msh")

# Save for running ABAQUS sim
write_abaqus_diffusion(sourcenodes,sinknodes)

# Visualise
gmsh.fltk.run()

# Finalise
gmsh.finalize()

# %%
"""
# Box with one central hole
"""

# %%
gmsh.initialize()

gmsh.model.add("t1")

gmsh.option.setNumber("Mesh.SaveGroupsOfElements", 1)
gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)

# Bounding box
boxsize = 1
gmsh.model.geo.addPoint(0, 0, 0, tag=1)
gmsh.model.geo.addPoint(boxsize, 0, 0, tag=2)
gmsh.model.geo.addPoint(boxsize, boxsize, 0, tag=3)
gmsh.model.geo.addPoint(0, boxsize, 0, tag=4)

gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(2, 3, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)

gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)

# Add 1 central circle to start with
r = 0.1
gmsh.model.geo.addPoint(boxsize/2, boxsize/2, 0, tag=5) # Centre
gmsh.model.geo.addPoint(boxsize/2, boxsize/2-r, 0, tag=6) # Rdial point south
gmsh.model.geo.addPoint(boxsize/2+r, boxsize/2, 0, tag=7) # Radial east
gmsh.model.geo.addPoint(boxsize/2, boxsize/2+r, 0, tag=8) # Radial north
gmsh.model.geo.addPoint(boxsize/2-r, boxsize/2, 0, tag=9) # Rdial point west
gmsh.model.geo.addCircleArc(6, 5, 7, 5)
gmsh.model.geo.addCircleArc(7, 5, 8, 6)
gmsh.model.geo.addCircleArc(8, 5, 9, 7)
gmsh.model.geo.addCircleArc(9, 5, 6, 8)
gmsh.model.geo.addCurveLoop([5,6,7,8], 2)

gmsh.model.geo.addPlaneSurface([1,2], 1) # Box with hole
gmsh.model.geo.addPlaneSurface([2], 2) # Circle

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(2, [1], name="layer0")
gmsh.model.addPhysicalGroup(2, [2], name="layer1")

# We can then generate a 2D mesh...
gmsh.model.mesh.generate(2)


sourcenodes = np.unique(gmsh.model.mesh.getElements(1,1)[-1][0])
sinknodes = np.unique(gmsh.model.mesh.getElements(1,3)[-1][0])
print(sourcenodes,sinknodes)

gmsh.write("gmshhole.msh")
write_abaqus_diffusion(sourcenodes,sinknodes)
gmsh.fltk.run()
gmsh.finalize()

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

gmsh.option.setNumber("Mesh.SaveGroupsOfElements", 1)
gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)
gmsh.option.setNumber("Geometry.OCCBooleanPreserveNumbering", 1)


# Bounding box
boxsize = 1
gmsh.model.occ.addRectangle(0,0,0,boxsize,boxsize,tag=1)

# Minimum spacing
eps = 1e-3

# Add a randomly centred circle
r = 0.2
ybuff = r + eps
while True:
    c = np.random.rand(2)*np.array([boxsize,boxsize-2*ybuff])
    c[1] += ybuff
    if np.abs(c[0]-r) > eps and np.abs(c[0]+r-boxsize) > eps:
        gmsh.model.occ.addDisk(c[0],c[1],0,r,r,tag=2)
        break

# Fragment overlapping shapes
frags, pc = gmsh.model.occ.fragment([(2, 1)], [(2, i) for i in range(2, 3)])
boxdimtag = pc[0][0]
boxtag = boxdimtag[1]

# Translate shapes in left periodic copy and refragment
eps = 1e-3
outs = gmsh.model.occ.getEntitiesInBoundingBox(-boxsize-eps,-eps,-eps,eps,boxsize+eps,eps,2)
gmsh.model.occ.translate(outs,boxsize,0,0)
# Check if bulk polymer entity has had label changed
frags, pc = gmsh.model.occ.fragment([(2, boxtag)], [(2, i) for i in range(2, len(frags)+1)])
boxdimtag = pc[0][0]
boxtag = boxdimtag[1]

# Translate shapes in right periodic copy and refragment
outs = gmsh.model.occ.getEntitiesInBoundingBox(boxsize-eps,-eps,-eps,2*boxsize+eps,boxsize+eps,eps,2)
gmsh.model.occ.translate(outs,-boxsize,0,0)
# Check if bulk polymer entity has had label changed
frags, pc = gmsh.model.occ.fragment([(2, boxtag)], [(2, i) for i in range(2, len(frags)+1)])
boxdimtag = pc[0][0]
boxtag = boxdimtag[1]

# Synchronise
gmsh.model.occ.synchronize()

#TODO Need to fix tag references to account for fragmenting
frags.remove(boxdimtag)
gmsh.model.addPhysicalGroup(2, [boxtag], name="layer0")
gmsh.model.addPhysicalGroup(2, [i[1] for i in frags], name="layer1")

# We can then generate a 2D mesh...
gmsh.model.mesh.generate(2)

gmsh.model.occ.synchronize()
gmsh.write("gmsh.msh")
sourceents = gmsh.model.occ.getEntitiesInBoundingBox(-eps,-eps,-eps,boxsize+eps,eps,eps,1)
sinkents = gmsh.model.occ.getEntitiesInBoundingBox(-eps,boxsize-eps,-eps,boxsize+eps,boxsize+eps,eps,1)
print(sourceents)
print(sinkents)
sourcenodes = np.unique(gmsh.model.mesh.getElements(1,sourceents[0][1])[-1][0])
sinknodes = np.unique(gmsh.model.mesh.getElements(1,sinkents[0][1])[-1][0])

gmsh.write("gmsh.msh")
write_abaqus_diffusion(sourcenodes,sinknodes)
#gmsh.fltk.run()
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

gmsh.option.setNumber("Mesh.SaveGroupsOfElements", 1)
gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)
gmsh.option.setNumber("Geometry.OCCBooleanPreserveNumbering", 1)

# Microstructure specifications
minspacefrac = 1e-3
vfrac = 0.5
r = 0.03
area0 = np.pi*r**2
print(area0)
nc = 100
area1 = nc*area0
print(area1)
area2 = area1/vfrac

# Bounding box
boxsize = np.sqrt(area2)
print(boxsize)
gmsh.model.occ.addRectangle(0,0,0,boxsize,boxsize,tag=1)

# Minimum spacing
eps = 1e-3*boxsize
print(eps*boxsize/area2)

# Add a randomly centred circle
ybuff = r + eps
cbuff = 2*r + eps
centers = np.empty((0,2))
print(cbuff)
# Loop over number of desired circles
for i in range(nc):

    # Continually try and insert circle within contstraints
    niter = 0
    reject = True
    while reject:

        # Randomly draw circle center
        c = np.random.rand(1,2)*np.array([[boxsize,boxsize-2*ybuff]])
        c[0,1] += ybuff
        
        # For first circle just check minimum distance from x-bounds
        if np.abs(c[0,0]-r) > eps and np.abs(c[0,0]+r-boxsize) > eps:
            if i == 0:
                reject = False
                        
            # Afterwards check against distance between previous circles
            else:
                reject = False
                for center in centers:
                    dist = np.linalg.norm(c[0]-center)
                    translator = np.array([boxsize,0.0])
                    ldist = np.linalg.norm(c[0]-translator-center)
                    rdist = np.linalg.norm(c[0]+translator-center)
                    mindist = np.min([dist,ldist,rdist])
                    if mindist < cbuff:
                        reject = True
        niter += 1
        if niter > 100000:
            raise Exception("Structure generation took too many iterations")

    # Add accepted circle
    print(i,c)
    centers = np.r_[centers,c]
    gmsh.model.occ.addDisk(c[0,0],c[0,1],0,r,r,tag=2+i)

centmod = centers.copy()
centmod *= 2
centmod[:,0] += 4.0
centmod[:,1] += 1.6
print(centmod)

# Fragment overlapping shapes
frags, pc = gmsh.model.occ.fragment([(2, 1)], [(2, i) for i in range(2, 2+nc)])
boxdimtag = pc[0][0]
boxtag = boxdimtag[1]

# Translate shapes in left periodic copy and refragment
eps = 1e-3
outs = gmsh.model.occ.getEntitiesInBoundingBox(-boxsize-eps,-eps,-eps,eps,boxsize+eps,eps,2)
gmsh.model.occ.translate(outs,boxsize,0,0)
# Check if bulk polymer entity has had label changed
ents = gmsh.model.occ.getEntities(2)
box = [boxdimtag]
ents.remove(boxdimtag)
frags, pc = gmsh.model.occ.fragment(box, ents)
boxdimtag = pc[0][0]
boxtag = boxdimtag[1]

# Translate shapes in right periodic copy and refragment
outs = gmsh.model.occ.getEntitiesInBoundingBox(boxsize-eps,-eps,-eps,2*boxsize+eps,boxsize+eps,eps,2)
gmsh.model.occ.translate(outs,-boxsize,0,0)
# Check if bulk polymer entity has had label changed
ents = gmsh.model.occ.getEntities(2)
box = [boxdimtag]
ents.remove(boxdimtag)
frags, pc = gmsh.model.occ.fragment(box, ents)
boxdimtag = pc[0][0]
boxtag = boxdimtag[1]

# Synchronise
gmsh.model.occ.synchronize()

#TODO Need to fix tag references to account for fragmenting
frags.remove(boxdimtag)
gmsh.model.addPhysicalGroup(2, [boxtag], name="layer0")
#gmsh.model.removeEntities([boxdimtag],True)
gmsh.model.addPhysicalGroup(2, [i[1] for i in frags], name="layer1")

# We can then generate a 2D mesh...
gmsh.option.setNumber("Mesh.MeshSizeMin",r)
gmsh.option.setNumber("Mesh.MeshSizeMax",r)
gmsh.model.mesh.generate(2)

gmsh.model.occ.synchronize()
sourceents = gmsh.model.occ.getEntitiesInBoundingBox(-eps,-eps,-eps,boxsize+eps,eps,eps,1)
sinkents = gmsh.model.occ.getEntitiesInBoundingBox(-eps,boxsize-eps,-eps,boxsize+eps,boxsize+eps,eps,1)
print(sourceents)
print(sinkents)
sourcenodes = np.unique(gmsh.model.mesh.getElements(1,sourceents[0][1])[-1][0])
sinknodes = np.unique(gmsh.model.mesh.getElements(1,sinkents[0][1])[-1][0])

gmsh.write("gmsh.msh")
write_abaqus_diffusion(sourcenodes,sinknodes)
gmsh.fltk.run()
gmsh.finalize()

# %%
r = 2
print(np.pi*r**2/(2*r)**2)
print(2*np.pi*r**2/(2*r)**2)
print(np.pi/4)

# %%
print(np.sqrt(r**2+r**2)-r)

# %%
