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
#import odbAccess

# Establish viewport
myViewport = session.Viewport(name='viewport',
    origin=(0.0, 0.0), width=679.98, height=459.45)
session.viewports['viewport'].makeCurrent()
session.viewports['viewport'].maximize()

# Permeation model
perm = mdb.Model(name='perm')

# Define rectangular sketch instance
sketch = perm.ConstrainedSketch(name='sketch',sheetSize=1.0)
sketch.rectangle(point1=(-0.5,0.5),point2=(0.5,-0.5))

# Define 2D shell part using sketch
part = perm.Part(dimensionality=TWO_D_PLANAR, name='part', type=DEFORMABLE_BODY)
part.BaseShell(sketch=sketch)

# Partition into two distinct layers
part.PartitionFaceByShortestPath(faces= \
    part.faces.findAt(((0.0,0.0,0.0), ), ), \
    point1=(0.0,-0.5,0.0), point2=(0.0,0.5,0.0))

# Define layer properties
layer1 = perm.Material(name='layer1')
layer1.Diffusivity(table=((1.0, 0.0), ))
layer1.diffusivity.setValues(law=FICK)
layer1.Solubility(table=((1.0, ), ))
layer2 = perm.Material(name='layer2')
layer2.Diffusivity(law=FICK, table=((0.1, 0.0), ))
layer2.Solubility(table=((1.1, ), ))

# Define and assign sections
perm.HomogeneousSolidSection(material='layer1', name='layer1',thickness=None)
perm.HomogeneousSolidSection(material='layer2', name='layer2',thickness=None)
setlayer1 = part.Set(faces=part.faces.findAt(((-0.25,0.0,0.0), ), ), name='layer1')
setlayer2 = part.Set(faces=part.faces.findAt(((0.25,0.0,0.0), ), ), name='layer2')
part.SectionAssignment(offset=0.0, \
    offsetField='', offsetType=MIDDLE_SURFACE, region=setlayer1, sectionName='layer1', \
    thicknessAssignment=FROM_SECTION)
part.SectionAssignment(offset=0.0, \
    offsetField='', offsetType=MIDDLE_SURFACE, region=setlayer2, sectionName='layer2', \
    thicknessAssignment=FROM_SECTION)

# Define assembly
#perm.rootAssembly.DatumCsysByDefault(CARTESIAN)
instance = perm.rootAssembly.Instance(dependent=OFF, name='instance', part=part)

# Set element types
#face1 = perm.rootAssembly.Set(faces=instance.faces.findAt(((-0.25,0.0,0.0), ), ), name='face1')
#face2 = perm.rootAssembly.Set(faces=instance.faces.findAt(((0.25,0.0,0.0), ), ), name='face2')
face1 = instance.faces.findAt(((-0.25,0.0,0.0), ), )
face2 = instance.faces.findAt(((0.25,0.0,0.0), ), )
perm.rootAssembly.setElementType(elemTypes=(ElemType( \
    elemCode=DC2D8, elemLibrary=STANDARD), ElemType(elemCode=DC2D6, \
    elemLibrary=STANDARD)), regions=(face1,face2, ))
    #elemLibrary=STANDARD)), regions=( \
    #part.faces.findAt(((-0.25,0.0,0.0), (0.25,0.0,0.0), ), ), ))

# Set mesh controls
perm.rootAssembly.setMeshControls(elemShape=QUAD, regions=face1, \
    technique=STRUCTURED)
perm.rootAssembly.setMeshControls(elemShape=QUAD, regions=face2, \
    technique=STRUCTURED)

# Seed mesh and generate
perm.rootAssembly.Set(edges= \
    instance.edges.findAt(((-0.25,0.5,0.0), ), ), name='top1')
perm.rootAssembly.Set(edges= \
    instance.edges.findAt(((0.25,0.5,0.0), ), ), name='top2')
perm.rootAssembly.seedEdgeByNumber(constraint=FIXED, \
    edges= instance.edges.findAt(((-0.25,0.5,0.0), ), ), \
    number=40)
perm.rootAssembly.seedEdgeByNumber(constraint=FIXED, \
    edges= instance.edges.findAt(((-0.25,-0.5,0.0), ), ), \
    number=40)
perm.rootAssembly.seedEdgeByNumber(constraint=FIXED, \
    edges= instance.edges.findAt(((0.25,0.5,0.0), ), ), \
    number=40)
perm.rootAssembly.seedEdgeByNumber(constraint=FIXED, \
    edges= instance.edges.findAt(((0.25,-0.5,0.0), ), ), \
    number=40)
perm.rootAssembly.seedEdgeByNumber(constraint=FIXED, edges= \
    instance.edges.findAt(((-0.5,0.0,0.0), ), ), number=1)
perm.rootAssembly.seedEdgeByNumber(constraint=FIXED, edges= \
    instance.edges.findAt(((0.5,0.0,0.0), ), ), number=1)
perm.rootAssembly.generateMesh(regions=(instance, ))

# Mass diffusion step
perm.MassDiffusionStep(dcmax=1.0, end=0.0, initialInc=0.00012, \
    maxInc=0.001, minInc=0.00012, name='diffusion', previous='Initial', \
    timePeriod=2.0, maxNumInc=100000)

# Concentration BCs
perm.rootAssembly.Set(edges= \
    instance.edges.findAt(((-0.5,0.0,0.0), ), ), name='left')
perm.rootAssembly.Set(edges= \
    instance.edges.findAt(((0.5,0.0,0.0), ), ), name='right')
perm.ConcentrationBC(amplitude=UNSET, createStepName=
    'diffusion', distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=
    1.0, name='source', region=perm.rootAssembly.sets['left'])
perm.ConcentrationBC(amplitude=UNSET, createStepName=
    'diffusion', distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=
    0.0, name='sink', region=perm.rootAssembly.sets['right'])

# Set absolute zero temperature
perm.setValues(absoluteZero=0.0)

# Output requests at specific times
perm.TimePoint(name='timepoints', points=((0.001, ), (0.05, 
    ), (0.2, ), (2.0, )))
perm.fieldOutputRequests['F-Output-1'].setValues(position=
    NODES, timePoint='timepoints', variables=('MFL', 'MFLT', 'CONC', 'NNC'))

# Create and submit job
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='perm', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name='sim', nodalOutputPrecision=SINGLE, 
    numCpus=1, numGPUs=0, numThreadsPerMpiProcess=1, queue=None, resultsFormat=
    ODB, scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, 
    waitMinutes=0)
#mdb.jobs['sim'].submit(datacheckJob=True)
#mdb.jobs['sim'].submit(continueJob=True)
mdb.jobs['sim'].submit()
mdb.jobs['sim'].waitForCompletion()
#mdb.jobs['sim'].submit()

# Open odb
o1 = session.openOdb(name='./sim.odb', \
    readOnly=False)
session.viewports['viewport'].setValues(displayedObject=o1)
odb = session.odbs['./sim.odb']

# Write field report to csv file
session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES)
session.writeFieldReport(fileName='check.csv', append=OFF,
    sortItem='Node Label', odb=odb, step=0, frame=0, outputPosition=NODAL,
    variable=(('CONC', ELEMENT_NODAL), ), stepFrame=ALL)

# Save model
mdb.saveAs(pathName='./model')
