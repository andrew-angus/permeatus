# %%
"""
# Initialisation
"""

# %%
import numpy as np
import matplotlib.pyplot as plt
import copy
import seaborn as sb
from cycler import cycler

# %%
%load_ext autoreload
%autoreload 2

# %%
# Matplotlib customisation with pgf output and nord colour scheme
nord0 = '#2E3440'
nord1 = '#3B4252'
nord2 = '#434C5E'
nord3 = '#4C566A'
nord4 = '#D8DEE9'
nord5 = '#E5E9F0'
nord6 = '#ECEFF4'
nord7 = '#8FBCBB'
nord8 = '#88C0D0'
nord9 = '#81A1C1'
nord10 = '#5E81AC'
nord11 = '#BF616A'
nord12 = '#D08770'
nord13 = '#EBCB8B'
nord14 = '#A3BE8C'
nord15 = '#B48EAD'
sb.set_theme(context='notebook',style='ticks')
plt.rcParams.update({
  "text.usetex": True,
  "text.latex.preamble": r'\usepackage{siunitx}',
  "font.family": 'serif',
  "figure.autolayout": True,
  "font.size": 11,
  "pgf.texsystem": "pdflatex",
  'pgf.rcfonts': False,
})
plt.rcParams['lines.linewidth'] = 1.5
w,h = plt.figaspect(1.618034)
textwidth = 5.50107
plt.rcParams['figure.figsize'] = [0.85*textwidth,0.618*0.85*textwidth]
plt.rc('text',usetex=True,color=nord0)
plt.rc('axes',edgecolor=nord0,labelcolor=nord0)
plt.rc('xtick',color=nord0)
plt.rc('ytick',color=nord0)
nord_cycler = (cycler(color=[nord10,nord11,nord7,nord3,nord15,nord12,nord13,nord14,nord0,nord8,nord9,nord4])\
               +cycler(linestyle=['-','--','-.',':','-','--',':','-.','-','--',':','-.']))
plt.rc('axes',prop_cycle=nord_cycler)

# %%
"""
# Permeatus 2-layer planar test
"""

# %%
from permeatus import *

# %%
perm = planar(layers=2,L=np.array([0.5,0.5]),D=np.array([1.0,0.1]),S=np.array([1.0,1.1]),\
              C0=1.0,C1=0.0,touts=[0.001,0.05,0.2,2.0],tstep=0.001,ncpu=1,N=[40,36])
#perm = planar(layers=2,L=np.array([0.5,0.5]),P=[1.0,0.11],\
#              p0=1.0,p1=0.0,touts=[0.001,0.05,0.2,2.0],tstep=0.001,ncpu=1,N=[40,36])

# %%
perm.submit_job()

# %%
perm.read_field()

# %%
perm.plot_1d()

# %%
xc, C, J = perm.steady_state('C',plot=True)

# %%
print(J)
print(C[1:-1])

# %%
xc, p, J = perm.steady_state('p',plot=True)

# %%
print(J)
print(p[1])

# %%
"""
# Permeatus n-layer planar
"""

# %%
#TODO use real polymer system
perm = planar(layers=3,L=np.array([0.4,0.3,0.3]),D=np.array([1.0,0.1,0.5]),S=np.array([1.0,1.1,0.7]),\
              C0=1.0,C1=0.0,touts=[0.001,0.05,0.2,2.0],tstep=0.001,ncpu=1,N=[40,36,32])

# %%
perm.submit_job()

# %%
perm.read_field()

# %%
perm.plot_1d('C')

# %%
perm.plot_1d('p')

# %%
xc, C, J = perm.steady_state('C',plot=True)

# %%
print(J)
print(C[1:-1])

# %%
xc, p, J = perm.steady_state('p',plot=True)

# %%
print(J)
print(p[1:-1])

# %%
perm.plot_1d('C',showplot=False)
scrap = perm.steady_state('C',plot=True,showplot=False)
plt.legend()
plt.show()

# %%
perm.plot_1d('p',showplot=False)
scrap = perm.steady_state('p',plot=True,showplot=False)
plt.legend()
plt.show()

# %%
