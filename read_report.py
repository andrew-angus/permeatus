import numpy as np
import csv

# Function which reads xy report from abaqus
def read_xy(fname='abaqus.rpt'):

  # Initialise data dict and labels logical
  data = {}
  labels = True

  # Open report file and go through lines
  with open(fname) as f:
    for line in f:
      # Filter lines which are not of interest (anything but data and labels)
      if '             ' in line and 'Legend' not in line:
        # Split resulting lines into list of relevant labels/data
        entries = line.split(' ')
        entries = [i for i in entries if not (i == '' or i == '\n')]

        # Check for NoValue and strip new line characters
        skip = [False for i in range(len(entries))]
        for i in range(len(entries)):
          entries[i] = entries[i].rstrip()
          if entries[i] == 'NoValue':
            entries[i] = 0.0
            #skip[i] = True

        # Store labels
        if labels:
          data = {i:np.empty(0) for i in entries}
          labels = False
        # Append data
        else:
          for i,label in enumerate(data):
            if not skip[i]:
              data[label] = np.append(data[label],float(entries[i]))

  return data

# Read field data output to csv file from abaqus
def read_field(fname='abaqus.csv'):

  # Initialise data dict and labels logical
  data = {}
  labels = True

  # Read into dictionary with csv module
  inc = -1; incc = -1
  with open(fname, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
      # Strip whitespace from keys
      row = {i.strip():j for i,j in zip(row.keys(),row.values())}

      # Only store keys of interest
      keys = ['Frame','X','Y','CONC']
      row = {i:row[i] for i in keys}

      # Split frame into increment and time and extract values
      splits = row['Frame'].split()
      increment = int(splits[1].strip(':'))
      time = float(splits[-1])

      # Check for new increment and initialise data dict
      if increment != inc:
        incc += 1
        data[incc] = {'t':time,'x':np.array([float(row['X'])]),\
            'y':np.array([float(row['Y'])]), 'C':np.array([float(row['CONC'])])}
        inc = increment
        print('inc',incc)

      # Else append data
      else:
        data[incc]['x'] = np.r_[data[incc]['x'],float(row['X'])]
        data[incc]['y'] = np.r_[data[incc]['y'],float(row['Y'])]
        data[incc]['C'] = np.r_[data[incc]['C'],float(row['CONC'])]

  # Store number of frames
  data['frames'] = incc+1

  return data
