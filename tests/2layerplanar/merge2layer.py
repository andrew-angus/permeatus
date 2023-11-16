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

# Variable parameters
r = 0.2
L = [0.5,0.5]
D = [1.0,0.1]
S = [1.0,1.1]
N = [40,40]
C0 = 1.0
C1 = 0.0
touts = [0.001,0.05,0.2,2.0]
tstep = 0.001
ncpu = 1

# Derivative parameters
totL = sum(L)
nlayer = len(L)

# Establish viewport
myViewport = session.Viewport(name='viewport',
    origin=(0.0, 0.0), width=679.98, height=459.45)
session.viewports['viewport'].makeCurrent()
session.viewports['viewport'].maximize()

# Permeation model
perm = mdb.Model(name='perm')

# Define rectangular sketch instance
ticker = 0.0
for i in range(nlayer):
  sketch = perm.ConstrainedSketch(name='sketch',sheetSize=max(2.1*totL,1.1*r))
  sketch.rectangle(point1=(0.0,ticker),point2=(r,ticker+L[i]))

  # Define 2D shell part using sketch
  perm.Part(dimensionality=TWO_D_PLANAR, name='part'+str(i), type=DEFORMABLE_BODY)
  perm.parts['part'+str(i)].BaseShell(sketch=sketch)
  ticker += L[i]

# Create assembly of part instances
instances = []
for i in range(nlayer):
  instance = perm.rootAssembly.Instance(dependent=OFF, name='part'+str(i)+str(i), 
      part=perm.parts['part'+str(i)])
  instances.append(instance)

# Merge all instances into new part
perm.rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
    instances=tuple(instances[i] for i in range(nlayer)), 
    keepIntersections=ON, name='part'+str(nlayer), originalInstances=SUPPRESS)
part = perm.parts['part'+str(nlayer)]

# Partition into distinct layers
#ticker = 0.0
#for i in range(nlayer-1):
#  ticker += L[i]
#  part.PartitionFaceByShortestPath(faces= \
#      part.faces.findAt(((ticker,0.5,0.0), ), ), \
#      point1=(ticker,0.0,0.0), point2=(ticker,1.0,0.0))

# Define and assign layer properties
ticker = 0.0
for i in range(nlayer):
  # Material properties
  perm.Material(name='layer'+str(i))
  perm.materials['layer'+str(i)].Diffusivity(table=((D[i], 0.0), ))
  perm.materials['layer'+str(i)].diffusivity.setValues(law=FICK)
  perm.materials['layer'+str(i)].Solubility(table=((S[i], ), ))

  # Material assignment to sections
  ticker += L[i]
  midlayer = ticker - L[i]/2
  perm.HomogeneousSolidSection(material='layer'+str(i), name='layer'+str(i),thickness=None)
  setlayer = part.Set(faces=part.faces.findAt(((r/2,midlayer,0.0), ), ), name='layer'+str(i))
  part.SectionAssignment(offset=0.0,offsetField='',offsetType=MIDDLE_SURFACE, \
      region=setlayer,sectionName='layer'+str(i),thicknessAssignment=FROM_SECTION)

# Redefine assembly
perm.rootAssembly.regenerate()
perm.rootAssembly.makeIndependent(instances=(
    perm.rootAssembly.instances['part'+str(nlayer)+'-1'], ))
instance = perm.rootAssembly.instances['part'+str(nlayer)+'-1']
#instance = perm.rootAssembly.Instance(dependent=OFF, name='instance', part=part)

# Mesh initialisation
face = []
ticker = 0.0
for i in range(nlayer):
  # Get face of layer
  ticker += L[i]
  midlayer = ticker - L[i]/2
  face.append(instance.faces.findAt(((r/2,midlayer,0.0), ), ))

  # Set mesh controls
  perm.rootAssembly.setMeshControls(elemShape=QUAD, regions=face[i], \
      technique=STRUCTURED)

# Set element types for all faces
perm.rootAssembly.setElementType(elemTypes=(ElemType( \
    elemCode=DC2D8, elemLibrary=STANDARD), ElemType(elemCode=DC2D6, \
    elemLibrary=STANDARD)), regions=tuple(face[i] for i in range(len(face))))

# Seed mesh edges and generate
# Left and right boundaries
perm.rootAssembly.seedEdgeByNumber(constraint=FIXED, edges= \
    instance.edges.findAt(((r/2,0.0,0.0), ), ), number=1)
perm.rootAssembly.seedEdgeByNumber(constraint=FIXED, edges= \
    instance.edges.findAt(((r/2,totL,0.0), ), ), number=1)

# Seed upper and lower edges for each layer
ticker = 0.0
maxele = 0.0
for i in range(nlayer):
  ticker += L[i]
  midlayer = ticker - L[i]/2
  perm.rootAssembly.seedEdgeByNumber(constraint=FIXED, \
      edges= instance.edges.findAt(((r,midlayer,0.0), ), ), \
      number=N[i])
  perm.rootAssembly.seedEdgeByNumber(constraint=FIXED, \
      edges= instance.edges.findAt(((0.0,midlayer,0.0), ), ), \
      number=N[i])

  # Track maximum element size in x-direction
  elesize = L[i]/N[i]
  if elesize > maxele:
    maxele = elesize

# Generate mesh
perm.rootAssembly.generateMesh(regions=(instance, ))

# Mass diffusion step
minInc = 1.01*maxele**2/(6*min(D))
maxInc = max(1.01*minInc,tstep) 
initInc = max(min(maxInc*0.99,tstep),minInc*1.001)
perm.MassDiffusionStep(dcmax=C0, end=0.0, initialInc=initInc, \
    maxInc=maxInc, minInc=minInc, name='diffusion', previous='Initial', \
    timePeriod=touts[-1], maxNumInc=int(touts[-1]/minInc)+1)

# Concentration BCs
perm.rootAssembly.Set(edges= \
    instance.edges.findAt(((r/2,0.0,0.0), ), ), name='left')
perm.rootAssembly.Set(edges= \
    instance.edges.findAt(((r/2,totL,0.0), ), ), name='right')
perm.ConcentrationBC(amplitude=UNSET, createStepName=
    'diffusion', distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=
    C0, name='source', region=perm.rootAssembly.sets['left'])
perm.ConcentrationBC(amplitude=UNSET, createStepName=
    'diffusion', distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=
    C1, name='sink', region=perm.rootAssembly.sets['right'])

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