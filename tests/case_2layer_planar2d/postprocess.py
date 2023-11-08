# Author: Andrew Angus

from odbAccess import openOdb
import sys
import pickle
import numpy as np
import inspect

# Open odb
odb = openOdb('./sim.odb')

step = odb.steps['diffusion']

#touts = np.array([0.001,0.05,0.2,2.0])
#conc_field = {}
#for i in range(len(touts)):
  #conc_field[i] = step.frames[i].fieldOutputs['CONC'].values
field = step.frames[-1].fieldOutputs['CONC']
for i in field.locations:
  print(i)
#for i in field.values:
#  print(i)
#  print(i.data)
  print('')
dp = field.values[0].data
print(dp)

#print(conc_field)

odb.close()

# Read frames for nodal concentrations
with open('test.pickle', 'wb') as handle:
    pickle.dump(dp, handle)

# Write field report to csv file
#session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES)
#session.writeFieldReport(fileName='check.csv', append=OFF,
#    sortItem='Node Label', odb=odb, step=0, frame=0, outputPosition=NODAL,
#    variable=(('CONC', ELEMENT_NODAL), ), stepFrame=ALL)

# Save model
#mdb.saveAs(pathName='./model')
