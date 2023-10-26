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

# Permeation model
perm = mdb.Model(name='perm')

# Define rectangular sketch instance
sketch = perm.ConstrainedSketch(name='sketch',sheetSize=1.0)
sketch.rectangle(point1=(-0.5,0.5),point2=(0.5,-0.5))

# Define 2D shell part using sketch
part = perm.Part(dimensionality=TWO_D_PLANAR, name='part', type=DEFORMABLE_BODY)
part.BaseShell(sketch=sketch)

# Partition into two distinct layers
print(part.faces)
part.PartitionFaceByShortestPath(faces= \
    part.faces.getSequenceFromMask(('[#1 ]', ), ), \
    point1=part.InterestingPoint( \
    part.edges[0], MIDDLE), point2= \
    part.InterestingPoint( \
    part.edges[2], MIDDLE))

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
setlayer1 = part.Set(faces=part.faces.getSequenceFromMask(('[#2 ]', ), ), name='layer1')
part.SectionAssignment(offset=0.0, \
    offsetField='', offsetType=MIDDLE_SURFACE, region=setlayer1, sectionName='layer1', \
    thicknessAssignment=FROM_SECTION)
setlayer2 = part.Set(faces=part.faces.getSequenceFromMask(('[#1 ]', ), ), name='layer2')
part.SectionAssignment(offset=0.0, \
    offsetField='', offsetType=MIDDLE_SURFACE, region=setlayer2, sectionName='layer2', \
    thicknessAssignment=FROM_SECTION)

# Define assembly
perm.rootAssembly.DatumCsysByDefault(CARTESIAN)
instance = perm.rootAssembly.Instance(dependent=OFF, name='instance', part=part)

# Mesh assembly instance
perm.rootAssembly.seedEdgeByBias(biasMethod=SINGLE, 
    constraint=FIXED, end2Edges=
    instance.edges.getSequenceFromMask(
    ('[#10 ]', ), ), number=40, ratio=5.0)
perm.rootAssembly.Set(edges=
    instance.edges.getSequenceFromMask(
    ('[#10 ]', ), ), name='top1')
perm.rootAssembly.seedEdgeByBias(biasMethod=SINGLE, 
    constraint=FIXED, end2Edges=
    instance.edges.getSequenceFromMask(
    ('[#8 ]', ), ), number=36, ratio=5.0)
perm.rootAssembly.Set(edges=
    instance.edges.getSequenceFromMask(
    ('[#8 ]', ), ), name='Edge Seeds-1')
perm.rootAssembly.seedEdgeByBias(biasMethod=SINGLE, 
    constraint=FIXED, end1Edges=
    instance.edges.getSequenceFromMask(
    ('[#2 ]', ), ), number=36, ratio=5.0)
perm.rootAssembly.seedEdgeByBias(biasMethod=SINGLE, 
    constraint=FIXED, end1Edges=
    instance.edges.getSequenceFromMask(
    ('[#40 ]', ), ), number=40, ratio=5.0)
perm.rootAssembly.seedEdgeByNumber(constraint=FIXED, edges=
    instance.edges.getSequenceFromMask(
    ('[#20 ]', ), ), number=1)
perm.rootAssembly.seedEdgeByNumber(constraint=FIXED, edges=
    instance.edges.getSequenceFromMask(
    ('[#4 ]', ), ), number=1)
perm.rootAssembly.generateMesh(regions=(
    instance, ))
perm.rootAssembly.setElementType(elemTypes=(ElemType(
    elemCode=DC2D8, elemLibrary=STANDARD), ElemType(elemCode=DC2D6, 
    elemLibrary=STANDARD)), regions=(
    instance.faces.getSequenceFromMask(
    ('[#3 ]', ), ), ))
perm.rootAssembly.deleteMesh(regions=
    instance.faces.getSequenceFromMask(
    ('[#3 ]', ), ))
perm.rootAssembly.setMeshControls(elemShape=QUAD, regions=
    instance.faces.getSequenceFromMask(
    ('[#3 ]', ), ), technique=STRUCTURED)
perm.rootAssembly.generateMesh(regions=(
    instance, ))

# Concentration BCs
perm.rootAssembly.Set(edges=
    instance.edges.getSequenceFromMask(
    ('[#20 ]', ), ), name='left')
perm.ConcentrationBC(amplitude=UNSET, createStepName=
    'diffusion', distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=
    1.0, name='source', region=perm.rootAssembly.sets['left'])
perm.rootAssembly.Set(edges=
    instance.edges.getSequenceFromMask(
    ('[#4 ]', ), ), name='right')
perm.ConcentrationBC(amplitude=UNSET, createStepName=
    'diffusion', distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=
    0.0, name='sink', region=perm.rootAssembly.sets['right'])

# Set absolute zero temperature
perm.setValues(absoluteZero=0.0)

# Mass diffusion step
perm.MassDiffusionStep(dcmax=1.0, end=0.0, initialInc=0.00012, 
    maxInc=0.001, minInc=0.00012, name='diffusion', previous='Initial', 
    timePeriod=2.0, maxNumInc=100000)

# Output requests at specific times
perm.TimePoint(name='timepoints', points=((0.001, ), (0.05, 
    ), (0.2, ), (2.0, )))
perm.fieldOutputRequests['fields'].setValues(position=
    NODES, timePoint='timepoints', variables=('MFL', 'MFLT', 'CONC', 'NNC'))

# Create and submit job
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='perm', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name='sim', nodalOutputPrecision=SINGLE, 
    numCpus=1, numGPUs=0, numThreadsPerMpiProcess=1, queue=None, resultsFormat=
    ODB, scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, 
    waitMinutes=0)
mdb.jobs['sim'].submit(consistencyChecking=OFF)
