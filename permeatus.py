import numpy as np
import csv
import matplotlib.pyplot as plt

# ABAQUS permeation class object
class permeatus:

  # Initialisation arguments
  def __init__(self,layers=None,L=None,D=None,S=None,P=None,\
               C0=None,C1=None,p0=None,p1=None,\
               N=None,tout=None,tstep=None,ncpu=None):

    #TODO Check validity of arguments 

    # Assign attributes
    self.layers = layers
    self.L = L
    self.D = D
    self.S = S
    self.P = P
    self.C0 = C0
    self.C1 = C1
    self.p0 = p0
    self.p1 = p1
    self.N = N
    self.tout = tout
    self.tstep = tstep
    self.ncpu = ncpu
    self.totL = np.sum(L)

    # Initialise calculated attributes
    self.J = None

    # Calculate P from DS, or vice versa
    if self.D is not None:
      self.P = self.D*self.S
    elif self.P is not None:
      self.D = self.P
      self.S = np.ones_like(self.D)

    # Calculate one BC from the other
    if self.p0 is not None:
      self.C0 = self.p0*self.S[0]
    elif self.C0 is not None:
      self.p0 = self.C0/self.S[0]
    if self.p1 is not None:
      self.C1 = self.p1*self.S[-1]
    elif self.C1 is not None:
      self.p1 = self.C1/self.S[-1]

  # Submit ABAQUS job
  #TODO options for resource control on job
  def submit_job(self):
    # Edit script template and run through abaqus interpreter
    pass

  # Read field data output to csv file from abaqus
  def read_field(self,fname='abaqus.csv'):

    # Initialise data dict and labels logical
    self.field = {}
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
          self.field[incc] = {'t':time,'x':np.array([float(row['X'])]),\
              'y':np.array([float(row['Y'])]), 'C':np.array([float(row['CONC'])])}
          inc = increment

        # Else append data
        else:
          self.field[incc]['x'] = np.r_[self.field[incc]['x'],float(row['X'])]
          self.field[incc]['y'] = np.r_[self.field[incc]['y'],float(row['Y'])]
          self.field[incc]['C'] = np.r_[self.field[incc]['C'],float(row['CONC'])]

    # Store number of frames
    self.field['frames'] = incc+1

  # Plot 1D solution
  def plot_1d(self,y='C'):

    # Loop through frames
    for frame in range(1,self.field['frames']):
      # Identify path along top of part
      ymax = np.max(self.field[frame]['y'])
      pathargs = np.ravel(np.argwhere(self.field[frame]['y'] > ymax-1e-10))
      x = self.field[frame]['x'][pathargs]
      C = self.field[frame]['C'][pathargs]

      # Sort by x-coordinate
      xsort = np.argsort(x)
      x = x[xsort]
      C = C[xsort]

      #TODO make sure interface points in correct order

      #TODO Plot either concentration or pressure
      plt.plot(x,C,label=f"{self.field[frame]['t']:0.3f} s")

    # Finalise plotting
    plt.legend()
    plt.xlabel(r'$x$ [$m$]')
    if y == 'C':
      plt.ylabel(r'$C$ [mol$m^{-3}$]')
    else:
      plt.ylabel(r'$p$ [Pa]')
    plt.show()

  # Function which reads xy report from abaqus
  def read_xy(self,fname='abaqus.rpt'):

    # Initialise data dict and labels logical
    self.xy = {}
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
            self.xy = {i:np.empty(0) for i in entries}
            labels = False
          # Append data
          else:
            for i,label in enumerate(self.xy):
              if not skip[i]:
                self.xy[label] = np.append(self.xy[label],float(entries[i]))

  # Linear algebra steady state solution
  def steady_state(self,y='C',plot=False,showplot=True):

    # Get linear coefficients relating pressure to molar flux
    k = self.P/self.L

    # Treat low layer number cases seperately
    # Evaluate pressure, concentration, and molar flux at relevant grid points
    if self.layers == 1:
      p = np.array([self.p0,self.p1])
      J = k[0]*(self.p0-self.p1)
      C = np.array([self.C0,self.C1])
      x = np.linspace(0,self.totL,2)
      xc = np.linspace(0,self.totL,2)
    elif self.layers == 2:
      p = np.zeros(3)
      p[0] = self.p0
      p[-1] = self.p1
      p[1] = k[0]*p[0]/(k[0]+k[1])
      J = k[0]*(p[0]-p[1])
      C = np.zeros(4)
      C[:2] = p[:2]*self.S[0]
      C[2:] = p[1:]*self.S[1]
      x = np.zeros_like(p)
      x[-1] = self.totL
      x[1] = self.L[0]
      xc = np.linspace(0,self.totL,4)
      xc[1:-1] = x[1]

    # n-layer solve
    else:
      # Populate internal pressure problem arrays
      # Incorporate boundary conditions
      b = np.zeros(self.layers-1)
      b[0] = self.p0
      b[-1] = self.p1
      # coefficient matrix
      a = -k[1:-1]*(np.eye(self.layers,k=1)+np.eye(self.layers,k=-1)) + \
          (k[:-1]+k[1:])*np.eye(self.layers)

      # Solve for internal pressures
      p_int = np.linalg.solve(a,b)

      # Calculate molar flux
      J = k[0]*(self.p0-p_int[0])

      # Construct full pressure array
      p = np.zeros(self.layers+1)
      p[0] = self.p0
      p[-1] = self.p1
      p[1:-1] = p_int

      # Get x points for pressures
      x = np.zeros_like(p)
      x[1:] += np.cumsum(self.L)

      # Calculate concentrations from pressures
      # x points
      xc = np.zeros(2*self.layers)
      xc[-1] = x[-1]
      xc[1:-2:2] = x[1:-1]
      xc[2:-1:2] = x[1:-1]
      # concentrations
      C = np.zeros_like(xc)
      C[:-1:2] = S*p[:-1]
      C[1::2] = S*p[1:]

    # Optionally plot
    if plot:
      if y == 'C':
        plt.plot(xc,C,label='steady-state')
        plt.ylabel(r'$C$ [mol$m^{-3}$]')
      else:
        plt.plot(x,p,label='steady-state')
        plt.ylabel(r'$p$ [Pa]')
      plt.xlabel(r'$x$ [$m$]')
      if showplot:
        plt.show()

    # Store and return results
    self.x = x
    self.p = p
    self.J = J
    self.xc = xc
    self.C = C

    if y == 'C':
      return xc, C, J
    else:
      return x, p, J



