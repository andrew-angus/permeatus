# Author: Andrew Angus

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from abaqus import *
from abaqusConstants import *
from math import *
import sys

# Variable parameters
inner_r = 0.5
dr = [0.5,0.5]
arc_extent = 0.25
D = [1.0,0.1]
S = [1.0,1.1]
N = [40,40]
C0 = 1.0
C1 = 0.0
touts = [0.001,0.05,0.2,2.0]
tstep = 0.001
ncpu = 1

# Derivative parameters
outer_r = sum(dr)+inner_r
nlayer = len(dr)

# Establish viewport
myViewport = session.Viewport(name='viewport',
    origin=(0.0, 0.0), width=679.98, height=459.45)
session.viewports['viewport'].makeCurrent()
session.viewports['viewport'].maximize()

# Permeation model
perm = mdb.Model(name='perm')

# Define sketch sheet
sketch = perm.ConstrainedSketch(name='sketch',sheetSize=2.1*outer_r)

# Get arc coordinates for inner and outer part limits
theta = arc_extent*pi
x1 = inner_r*cos(theta)
y1 = inner_r*sin(theta)
x2 = inner_r*cos(-theta)
y2 = inner_r*sin(-theta)
x3 = outer_r*cos(theta)
y3 = outer_r*sin(theta)
x4 = outer_r*cos(-theta)
y4 = outer_r*sin(-theta)

# Draw part sketch
sketch.ArcByCenterEnds(center=(0.0, 0.0),direction=CLOCKWISE, \
    point1=(x1, y1), point2=(x2, y2))
sketch.ArcByCenterEnds(center=(0.0, 0.0),direction=CLOCKWISE, \
    point1=(x3, y3), point2=(x4, y4))
sketch.Line(point1=(x1, y1), point2=(x3, y3))
sketch.Line(point1=(x2, y2), point2=(x4, y4))

# Define 3D part by extruding sketch and rounding off corners till circular
part = perm.Part(dimensionality=TWO_D_PLANAR, name='part', type=DEFORMABLE_BODY)
part.BaseShell(sketch=sketch)

mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-1', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-1'].BaseSolidExtrude(depth=0.2, sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
mdb.models['Model-1'].parts['Part-1'].Round(edgeList=(
    mdb.models['Model-1'].parts['Part-1'].edges[4], 
    mdb.models['Model-1'].parts['Part-1'].edges[6], 
    mdb.models['Model-1'].parts['Part-1'].edges[10], 
    mdb.models['Model-1'].parts['Part-1'].edges[11]), radius=0.1)


# Partition into distinct layers
ticker = inner_r
for i in range(nlayer-1):
  ticker += dr[i]
  x1 = ticker*cos(theta)
  y1 = ticker*sin(theta)
  x2 = ticker*cos(-theta)
  y2 = ticker*sin(-theta)
  part.PartitionFaceByCurvedPathEdgePoints(\
      edge1=part.edges.findAt(((x1,y1,0.0), ), )[0], \
      edge2=part.edges.findAt(((x2,y2,0.0), ), )[0], \
      face=part.faces.findAt(((ticker,0.0,0.0), ), )[0], \
      point1=(x1,y1,0.0),point2=(x2,y2,0.0))

# Define and assign layer properties
ticker = inner_r
for i in range(nlayer):
  # Material properties
  perm.Material(name='layer'+str(i))
  perm.materials['layer'+str(i)].Diffusivity(table=((D[i], 0.0), ))
  perm.materials['layer'+str(i)].diffusivity.setValues(law=FICK)
  perm.materials['layer'+str(i)].Solubility(table=((S[i], ), ))

  # Material assignment to sections
  ticker += dr[i]
  midlayer = ticker - dr[i]/2
  perm.HomogeneousSolidSection(material='layer'+str(i), name='layer'+str(i),thickness=None)
  setlayer = part.Set(faces=part.faces.findAt(((midlayer,0.0,0.0), ), ), name='layer'+str(i))
  part.SectionAssignment(offset=0.0,offsetField='',offsetType=MIDDLE_SURFACE, \
      region=setlayer,sectionName='layer'+str(i),thicknessAssignment=FROM_SECTION)

# Define assembly
instance = perm.rootAssembly.Instance(dependent=OFF, name='instance', part=part)

# Mesh initialisation
face = []
ticker = inner_r
for i in range(nlayer):
  # Get face of layer
  ticker += dr[i]
  midlayer = ticker - dr[i]/2
  face.append(instance.faces.findAt(((midlayer,0.0,0.0), ), ))

  # Set mesh controls
  perm.rootAssembly.setMeshControls(elemShape=TRI, regions=face[i], \
      technique=STRUCTURED)

# Set element types for all faces
perm.rootAssembly.setElementType(elemTypes=(ElemType( \
    elemCode=DC2D8, elemLibrary=STANDARD), ElemType(elemCode=DC2D6, \
    elemLibrary=STANDARD)), regions=tuple(face[i] for i in range(len(face))))

mdb.saveAs(pathName='./model')
#sys.exit()

# Seed mesh edges and generate
# Inner and outer boundaries
perm.rootAssembly.seedEdgeByNumber(constraint=FIXED, edges= \
    instance.edges.findAt(((inner_r,0.0,0.0), ), ), number=max(1,int(theta/(pi*9/360))))
perm.rootAssembly.seedEdgeByNumber(constraint=FIXED, edges= \
    instance.edges.findAt(((outer_r,0.0,0.0), ), ), number=max(1,int(theta/(pi*9/360))))

# Seed upper and lower edges for each layer
ticker = inner_r
maxele = 0.0
for i in range(nlayer):
  ticker += dr[i]
  midlayer = ticker - dr[i]/2
  x1 = midlayer*cos(theta)
  y1 = midlayer*sin(theta)
  x2 = midlayer*cos(-theta)
  y2 = midlayer*sin(-theta)
  perm.rootAssembly.seedEdgeByNumber(constraint=FIXED, \
      edges= instance.edges.findAt(((x1,y1,0.0), ), ), \
      number=N[i])
  perm.rootAssembly.seedEdgeByNumber(constraint=FIXED, \
      edges= instance.edges.findAt(((x2,y2,0.0), ), ), \
      number=N[i])

  # Track maximum element size in r-direction
  elesize = dr[i]/N[i]
  if elesize > maxele:
    maxele = elesize

# Generate mesh
perm.rootAssembly.generateMesh(regions=(instance, ))

# Mass diffusion step
minInc = 1.05*maxele**2/(6*min(D))
maxInc = max(minInc,tstep) 
initInc = min(maxInc,tstep)
perm.MassDiffusionStep(dcmax=C0, end=0.0, initialInc=initInc, \
    maxInc=maxInc, minInc=minInc, name='diffusion', previous='Initial', \
    timePeriod=touts[-1], maxNumInc=int(touts[-1]/minInc)+1)


# Concentration BCs
perm.rootAssembly.Set(edges= \
    instance.edges.findAt(((inner_r,0.0,0.0), ), ), name='inner')
perm.rootAssembly.Set(edges= \
    instance.edges.findAt(((outer_r,0.0,0.0), ), ), name='outer')
perm.ConcentrationBC(amplitude=UNSET, createStepName=
    'diffusion', distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=
    C0, name='source', region=perm.rootAssembly.sets['inner'])
perm.ConcentrationBC(amplitude=UNSET, createStepName=
    'diffusion', distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=
    C1, name='sink', region=perm.rootAssembly.sets['outer'])

# Set absolute zero temperature
perm.setValues(absoluteZero=0.0)

# Output requests at specific times
tpoints = tuple(tuple(touts[i] for j in range(1)) for i in range(len(touts)))
perm.TimePoint(name='timepoints', points=tpoints)
perm.fieldOutputRequests['F-Output-1'].setValues(position=
    NODES, timePoint='timepoints', variables=('MFL', 'MFLT', 'CONC', 'NNC'))

# Create and submit job
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='perm', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name='sim', nodalOutputPrecision=SINGLE, 
    numCpus=ncpu, numGPUs=0, numThreadsPerMpiProcess=1, queue=None, resultsFormat=
    ODB, scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, 
    waitMinutes=0)
mdb.saveAs(pathName='./model')
mdb.jobs['sim'].submit()
mdb.jobs['sim'].waitForCompletion()

# Open odb
o1 = session.openOdb(name='./sim.odb', \
    readOnly=False)
session.viewports['viewport'].setValues(displayedObject=o1)
odb = session.odbs['./sim.odb']

# Write concentration field report to csv file
session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES)
session.writeFieldReport(fileName='CONC.csv', append=OFF,
    sortItem='Node Label', odb=odb, step=0, frame=0, outputPosition=NODAL,
    variable=(('CONC', ELEMENT_NODAL), ), stepFrame=ALL)

# Save model
mdb.saveAs(pathName='./model')
