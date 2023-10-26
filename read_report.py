import numpy as np

def read_report(fname='abaqus.rpt'):

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



