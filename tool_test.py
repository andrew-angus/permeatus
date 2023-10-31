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
# Read xy report
"""

# %%
# Read 1 layer abaqus report
from read_report import *
data = read_xy('1layer.rpt')
#print(data)

# %%
# 1 layer plat
touts = [0.001,0.05,0.2,2.0]
for i in range(4):
    plt.plot(data['X']+0.5,data[f'tout{i+1}'],label=f'{touts[i]:0.3f} s')
plt.xlabel('$x$ [$m$]')
plt.ylabel(r'$\varphi$ [mol/$m^3$]')
plt.legend()
plt.tight_layout()
plt.savefig('1layer_abaqus.pgf',bbox_inches='tight')
plt.show()

# %%
# 2 layers
data = read_xy('2layer.rpt')
for i in range(4):
    datx = data['X'][np.nonzero(data[f'tout{i+1}'])]+0.5
    sortargs = np.argsort(datx)
    datx = datx[sortargs]
    daty = data[f'tout{i+1}'][np.nonzero(data[f'tout{i+1}'])]
    daty = daty[sortargs]
    datx = np.r_[datx,1.0]
    daty = np.r_[daty,0.0]
    plt.plot(datx,daty,label=f'{touts[i]:0.3f} s')
plt.xlabel('$x$ [$m$]')
plt.ylabel(r'$\varphi$ [mol/$m^3$]')
plt.legend()
plt.tight_layout()
plt.savefig('2layer_abaqus.pgf',bbox_inches='tight')
plt.show()

# %%
"""
# Permeatus 2-layer planar test
"""

# %%
from permeatus import *

# %%
perm = permeatus(layers=2,L=np.array([0.5,0.5]),D=np.array([1.0,0.1]),S=np.array([1.0,1.1]),C0=1.0,C1=0.0)

# %%
perm.read_field('./case_2layer_planar2d/check.csv')

# %%
perm.plot_1d()

# %%
a = np.random.rand(10)
print(a)
print(a[1:-2:2])

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
