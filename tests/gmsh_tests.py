# %%
"""
# ABAQUS 2-layer
"""

# %%
import gmsh
import sys
import fileinput as fi

gmsh.initialize()

gmsh.model.add("t1")

gmsh.option.setNumber("Mesh.SaveGroupsOfElements", 1)
gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)

width = 1.0/80
gmsh.model.geo.addPoint(0, 0, 0, tag=1,meshSize=width)
gmsh.model.geo.addPoint(width*10, 0, 0, tag=2)#,meshSize=width)
gmsh.model.geo.addPoint(width*10, .5, 0, tag=3)#,meshSize=width)
gmsh.model.geo.addPoint(0, .5, 0, tag=4)#,meshSize=width)
gmsh.model.geo.addPoint(width*10, 1.0, 0, tag=5)#,meshSize=width)
gmsh.model.geo.addPoint(.0, 1.0, 0, tag=6)#,meshSize=width)

gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(2, 3, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)

gmsh.model.geo.addLine(3, 5, 5)
gmsh.model.geo.addLine(5, 6, 6)
gmsh.model.geo.addLine(6, 4, 7)

# Need to both be anti-clockwise, or both clockwise
gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
gmsh.model.geo.addCurveLoop([-3, 5, 6, 7], 2)

gmsh.model.geo.addPlaneSurface([1], 1)
gmsh.model.geo.addPlaneSurface([2], 2)

gmsh.model.geo.mesh.setTransfiniteCurve(1, 2, "Progression", 1.0)
gmsh.model.geo.mesh.setTransfiniteCurve(6, 2)
gmsh.model.geo.mesh.setTransfiniteCurve(3, 2)
gmsh.model.geo.mesh.setTransfiniteCurve(2, 40)
gmsh.model.geo.mesh.setTransfiniteCurve(4, 40)
gmsh.model.geo.mesh.setTransfiniteCurve(5, 40)
gmsh.model.geo.mesh.setTransfiniteCurve(7, 40)
gmsh.model.geo.mesh.setTransfiniteSurface(1)
gmsh.model.geo.mesh.setTransfiniteSurface(2)

gmsh.model.geo.synchronize()

# gmsh.model.addPhysicalGroup(1, [1, 2, 4], 5)
#gmsh.model.addPhysicalGroup(1, [1], name="source")
#gmsh.model.addPhysicalGroup(1, [6], name="sink")
gmsh.model.addPhysicalGroup(2, [1], name="layer0")
gmsh.model.addPhysicalGroup(2, [2], name="layer1")

# We can then generate a 2D mesh...
#print(gmsh.option.getNumber("Mesh.MeshSizeMin"))
#gmsh.option.setNumber("Mesh.MeshSizeMin",width)
#print(gmsh.option.getNumber("Mesh.MeshSizeMin"))
#print(gmsh.option.getNumber("Mesh.MeshSizeMax"))
#gmsh.option.setNumber("Mesh.MeshSizeMax",width)
#print(gmsh.option.getNumber("Mesh.MeshSizeMax"))
#gmsh.model.mesh.setSizeAtParametricPoints(2,1,[0.0,0.0,width,0.0],[0.01])
#gmsh.model.mesh.setSizeAtParametricPoints(2,2,[0.width.0,width,1.0],[0.01])
#gmsh.model.mesh.setAlgorithm(2,1,8)
#gmsh.model.mesh.setAlgorithm(2,2,8)
#gmsh.model.mesh.setSizeAtParametricPoints(2,1,[width,0.0],[0.01])
#gmsh.model.mesh.embed(0,[7],2,1)
#gmsh.model.mesh.embed(1,[1],2,1)
#gmsh.model.mesh.embed(0,[7,8],2,1)
#gmsh.model.geo.addLine(7, 8, 8)
#gmsh.model.geo.synchronize()
#gmsh.model.mesh.embed(1,[8],2,1)
gmsh.model.mesh.generate(2)
gmsh.model.mesh.recombine()

gmsh.fltk.run()

sourcenodes = gmsh.model.mesh.getElements(1,1)[-1][0]
sinknodes = gmsh.model.mesh.getElements(1,6)[-1][0]
#print(gmsh.model.mesh.getElements())

gmsh.write("gmsh.inp")

# Replace element types with diffusion
with fi.input("gmsh.inp",inplace=True) as f:
    for line in f:
        print(line.replace("CPS4","DC2D4"), end='')

# Append BC's, Diffusion step details, Material properties and section assignment
with open("gmsh.inp","a") as f:
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
with open("gmsh.inp","r") as f:
    print(f.read())

gmsh.finalize()

# %%
"""
# RVE
"""

# %%
import gmsh
import sys
import fileinput as fi
import numpy as np

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

# Random circle insertion
vfrac = 0.5

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


# Mesh
"""
gmsh.model.geo.mesh.setTransfiniteCurve(1, 2, "Progression", 1.0)
gmsh.model.geo.mesh.setTransfiniteCurve(6, 2)
gmsh.model.geo.mesh.setTransfiniteCurve(3, 2)
gmsh.model.geo.mesh.setTransfiniteCurve(2, 40)
gmsh.model.geo.mesh.setTransfiniteCurve(4, 40)
gmsh.model.geo.mesh.setTransfiniteCurve(5, 40)
gmsh.model.geo.mesh.setTransfiniteCurve(7, 40)
gmsh.model.geo.mesh.setTransfiniteSurface(1)
gmsh.model.geo.mesh.setTransfiniteSurface(2)
"""

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(2, [1], name="organic")
gmsh.model.addPhysicalGroup(2, [2], name="inorganic")

# We can then generate a 2D mesh...
gmsh.model.mesh.generate(2)
#gmsh.model.mesh.recombine()

gmsh.fltk.run()

sourcenodes = np.unique(gmsh.model.mesh.getElements(1,1)[-1][0])
sinknodes = np.unique(gmsh.model.mesh.getElements(1,3)[-1][0])
#print(gmsh.model.mesh.getElements())
print(sourcenodes,sinknodes)

gmsh.write("gmsh.inp")

# Replace element types with diffusion
with fi.input("gmsh.inp",inplace=True) as f:
    for line in f:
        print(line.replace("CPS","DC2D"), end='')

# Append BC's, Diffusion step details, Material properties and section assignment
with open("gmsh.inp","a") as f:
    # Define boundary node sets
    f.write('*NSET, NSET=source \n')
    for i,j in enumerate(sourcenodes):
        if (i+1) % 10 == 0 or i + 1 == len(sourcenodes):
            f.write(f'{j},\n')
        else:
            f.write(f'{j}, ')
    f.write('*NSET, NSET=sink \n')
    for i,j in enumerate(sinknodes):
        if (i+1) % 10 == 0 or i + 1 == len(sinknodes):
            f.write(f'{j},\n')
        else:
            f.write(f'{j}, ')

    # Define materials
    f.write(f'*Material, name=organic\n')
    f.write(f'*Diffusivity, law=FICK\n')
    f.write(f'1.,0.\n')
    f.write(f'*Solubility\n')
    f.write(f'1.,\n')
    f.write(f'*Material, name=inorganic\n')
    f.write(f'*Diffusivity, law=FICK\n')
    f.write(f'0.1,0.\n')
    f.write(f'*Solubility\n')
    f.write(f'1.1,\n')

    # Assign material sections
    f.write(f'*Solid Section, elset=organic, material=organic\n')
    f.write(f'*Solid Section, elset=inorganic, material=inorganic\n')

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
with open("gmsh.inp","r") as f:
    print(f.read())

gmsh.finalize()

# %%
