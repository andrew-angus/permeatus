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
from permeatus import *

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
  #"text.latex.preamble": r'\usepackage{siunitx}',
  "font.family": 'serif',
  "figure.autolayout": True,
  "font.size": 11,
  #'mathtext.fontset': 'cm',
  "pgf.texsystem": "pdflatex",
  'pgf.rcfonts': False,
})
plt.rcParams['lines.linewidth'] = 1.5
w,h = plt.figaspect(1.618034)
textwidth = 5.50107
#textwidth = 5.5129
plt.rcParams['figure.figsize'] = [0.85*textwidth,0.618*0.85*textwidth]
#plt.rcParams['figure.figsize'] = [0.75*0.49*textwidth,0.75*0.49*textwidth]
plt.rc('text',usetex=True,color=nord0)
plt.rc('axes',edgecolor=nord0,labelcolor=nord0)
plt.rc('xtick',color=nord0)
plt.rc('ytick',color=nord0)
nord_cycler = (cycler(color=[nord10,nord11,nord7,nord3,nord15,nord12,nord13,nord14,nord0,nord8,nord9,nord4])\
               +cycler(linestyle=['-','--','-.',':','-','--','-.',':','-','--','-.',':']))
plt.rc('axes',prop_cycle=nord_cycler)

# %%
"""
# Permeatus 2-layer planar test
"""

# %%
perm = planar(layers=2,r=0.1,L=np.array([0.5,0.5]),D=np.array([1.0,0.1]),S=np.array([1.0,1.1]),\
              C0=1.0,C1=0.0,touts=[0.001,0.05,0.2,2.0],tstep=0.001,ncpu=1,N=[40,36])
#perm = planar(layers=2,L=np.array([0.5,0.5]),P=[1.0,0.11],\
#              p0=1.0,p1=0.0,touts=[0.001,0.05,0.2,2.0],tstep=0.001,ncpu=1,N=[40,36])

# %%
perm.submit_job()

# %%
perm.read_field('C')

# %%
perm.read_field('J')

# %%
print(perm.field)

# %%
print(np.mean(perm.field['J'][-1]['data'][:,1]))
print(np.std(perm.field['J'][-1]['data'][:,1]))
print(perm.field['J'][-1]['data'])

# %%
perm.plot_1d(showplot=False)
plt.tight_layout()
plt.savefig('abaqus2layer.pdf',bbox_inches='tight')
plt.show()

# %%
print(perm.frames)

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
perm = planar(layers=3,r=0.1,L=np.array([0.4,0.3,0.3]),D=np.array([1.0,0.1,0.5]),S=np.array([1.0,1.1,0.7]),\
              C0=1.0,C1=0.0,touts=[0.001,0.05,0.2,2.0],tstep=0.001,ncpu=1,N=[40,36,32])
#perm = planar(layers=3,L=np.array([8.5e-3,5.5e-3,6e-3]),P=np.array([1.0,0.1,0.5]),\
#             p0=35e6,p1=0.0,touts=[0.001,0.05,0.2,2.0],tstep=0.001,ncpu=1,N=[40,36,32])

# %%
perm.submit_job()

# %%
perm.read_field('C')

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
x, p, J = perm.steady_state('C',plot=True,showplot=False)
plt.legend()
plt.tight_layout()
plt.savefig('abaqus3layer.pgf',bbox_inches='tight')
plt.show()

# %%
perm.plot_1d('p',showplot=False)
scrap = perm.steady_state('p',plot=True,showplot=False)
plt.legend()
plt.tight_layout()
plt.savefig('abaqus3layerp.pdf',bbox_inches='tight')
plt.show()

# %%
"""
# Nielsen Model
"""

# %%
#TODO use real polymer system
permniels = planar(layers=3,r=0.1,L=np.array([0.4,0.3,0.3]),Dc=np.array([1.0,0.1,0.5]),Sc=np.array([1.0,1.1,0.7]),\
              C0=1.0,C1=0.0,touts=[0.001,0.05,0.2,2.0],tstep=0.001,ncpu=1,N=[40,36,32],Vd_frac=np.array([0.0,0.6,0.0]),AR=np.array([1.0,2.0,1.0]),model='Nielsen')

# %%
print(permniels.D)
print(permniels.S)
print(permniels.P)

# %%
permniels.submit_job()

# %%
permniels.read_field()

# %%
permniels.plot_1d('C',showplot=False)
x, p, J = permniels.steady_state('C',plot=True,showplot=False)
plt.legend()
plt.show()

# %%
print(J)
print(C[1:-1])

# %%
permniels.plot_1d('p',showplot=False,timemask=[False,False,False,False,True],plotlabels=[None,None,None,None,r'Central layer $\phi_{d} = 0.6$'])
x, p, J = perm.steady_state('p',plot=True,showplot=False,plotlabel='Homogeneous')
plt.legend()
plt.tight_layout()
plt.savefig('abaqusnielsen.pdf',bbox_inches='tight')
plt.show()

# %%
print(J)
print(C[1:-1])

# %%
print(0.082/0.251)

# %%
perm.plot_1d('p',showplot=False)
permniels.plot_1d('p',showplot=False)
plt.legend()
#plt.tight_layout()
#plt.savefig('abaqusnielsen.pdf',bbox_inches='tight')
plt.show()

# %%
print(J)
print(p[1:-1])

# %%
"""
## Real System
"""

# %%
"""
## Average coefficients
"""

# %%
perm.get_avg_coeffs()
print(perm.Davg,perm.Savg,perm.Pavg)

# %%
# Check above with single layer average sim
permavg = planar(layers=1,r=0.1,L=np.array([1.0]),D=np.array([perm.Davg]),S=np.array([perm.Savg]),\
              C0=1.0,C1=0.0,touts=[0.001,0.05,0.2,2.5],tstep=0.001,ncpu=1,N=[80])

# %%
permavg.submit_job()

# %%
permavg.read_field()

# %%
permavg.plot_1d('C',showplot=False)
x, p, J = permavg.steady_state('C',plot=True,showplot=False)
plt.legend()
plt.show()

# %%
print(J)

# %%
permavg.plot_1d('p',showplot=False)
xc, C, J = permavg.steady_state('p',plot=True,showplot=False)
plt.legend()
plt.tight_layout()
#plt.savefig('abaqusnielsen.pgf',bbox_inches='tight')
plt.show()

# %%
print(J)

# %%
"""
# Effective property bounds
"""

# %%
permbnd = planar(layers=3,L=np.array([0.4,0.3,0.3]),Dc=np.array([1.0,0.1,0.5]),Sc=np.array([1.0,1.1,0.7]),\
              C0=1.0,C1=0.0,touts=[0.001,0.05,0.2,2.0],tstep=0.001,ncpu=1,N=[40,36,32],Vd_frac=np.array([0.0,0.5,0.0]),AR=np.array([1.0,2.0,1.0]),\
              Dd=np.array([0.0,0.1,0.0]),Sd=np.array([0.0,1.0,0.0]),model=None)

# %%
# Get bounds by Wiener
permbnd.get_eff_bnds(method='Wiener',set_eff='false')
print(permbnd.D_lower,permbnd.D_upper)
print(permbnd.S_lower,permbnd.S_upper)
print(permbnd.P_lower,permbnd.P_upper)

# %%
# Compare two layer averaging with bounds method
permbnd = planar(layers=2,L=np.array([0.4,0.6]),Dc=np.array([1.0,0.1]),Sc=np.array([1.0,1.1]),\
              C0=1.0,C1=0.0,touts=[0.001,0.05,0.2,2.0],tstep=0.001,ncpu=1,N=[40,36],Vd_frac=np.array([0.4,0.6]),AR=np.array([1.0,2.0]),\
              Dd=np.array([0.1,1.0]),Sd=np.array([1.0,1.1]),model=None)

# %%
# Get bounds by Wiener
permbnd.get_eff_bnds(method='Wiener',set_eff='false')
print(permbnd.D_lower,permbnd.D_upper)
print(permbnd.S_lower,permbnd.S_upper)
print(permbnd.P_lower,permbnd.P_upper)

# %%
xp, p, J = permbnd.steady_state('p',plot=True)
permbnd.get_avg_coeffs()
print(permbnd.Davg)
print(permbnd.Savg)
print(permbnd.Pavg)

# %%
# Compare two layer averaging with bounds method
permbnd = planar(layers=2,L=np.array([0.4,0.6]),D=np.array([1.0,0.1]),S=np.array([1.0,1.1]),\
              C0=1.0,C1=0.0,touts=[0.001,0.05,0.2,2.0],tstep=0.001,ncpu=1,N=[40,36])

# %%
xp, p, J = permbnd.steady_state('p',plot=True)
permbnd.get_avg_coeffs()
print(permbnd.Davg)
print(permbnd.Savg)
print(permbnd.Pavg)
print(J)

# %%
Du = 0.4*1+0.6*0.1
Su = 0.4*1+0.6*1.1
Dl = 1/(0.4/1+0.6/0.1)
Sl = 1/(0.4/1+0.6/1.1)

# %%
print(Du,Su)
print(Dl,Sl)

# %%
print(Du*Su)
print(Dl*Sl)

# %%
# Check above with single layer average sim
permavg = planar(layers=1,r=0.1,L=np.array([1.0]),D=np.array([Du]),S=np.array([Su]),\
              C0=1.0,C1=0.0,touts=[0.001,0.05,0.2,2.5],tstep=0.001,ncpu=1,N=[80])

# %%
xp, p, J = permavg.steady_state('p',plot=True)
permavg.get_avg_coeffs()
print(permavg.Davg)
print(permavg.Savg)
print(permavg.Pavg)
print(J)

# %%
# Check above with single layer average sim
permavg = planar(layers=1,r=0.1,L=np.array([1.0]),D=np.array([Dl]),S=np.array([Sl]),\
              C0=1.0,C1=0.0,touts=[0.001,0.05,0.2,2.5],tstep=0.001,ncpu=1,N=[80])

# %%
xp, p, J = permavg.steady_state('C',plot=True)
permavg.get_avg_coeffs()
print(permavg.Davg)
print(permavg.Savg)
print(permavg.Pavg)
print(J)

# %%
"""
## Reproduce paper figure
"""

# %%
# Reproduce ebermannAnalytical2022 figure of lower and upper bounds vs volume composite volume fraction
wienerl = np.empty(0)
wieneru = np.empty(0)
HASl = np.empty(0)
HASu = np.empty(0)
points = 100
Vd_frac = np.linspace(0,1,points)
for i in range(points):
    permbnd = planar(layers=1,Pc=np.array([1.4e-17]),Pd=np.array([1.49e-18]),Vd_frac=np.array([Vd_frac[i]]))
    permbnd.get_eff_bnds('Wiener')
    wienerl = np.r_[wienerl,permbnd.P_lower/permbnd.Pc]
    wieneru = np.r_[wieneru,permbnd.P_upper/permbnd.Pc]
    permbnd.get_eff_bnds('HS')
    HASl = np.r_[HASl,permbnd.P_lower/permbnd.Pc]
    HASu = np.r_[HASu,permbnd.P_upper/permbnd.Pc]

# %%
plt.plot(Vd_frac,wienerl,label='Wiener lower')
plt.plot(Vd_frac,wieneru,label='Wiener upper')
plt.plot(Vd_frac,HASl,label='HS lower')
plt.plot(Vd_frac,HASu,label='HS upper')
plt.legend()
plt.grid()
plt.xlabel(r'$\phi_d$')
#plt.ylabel(r'$\frac{P_{eff}}{P_c}$',rotation=0,labelpad=15)
plt.ylabel(r'$\overline{P}/P_c$')#,rotation=0,labelpad=15)
plt.tight_layout()
plt.savefig('bounds.pdf',bbox_inches='tight')
plt.show()

# %%
print(permbnd.Pd/permbnd.Pc)
print(permbnd.Pc/permbnd.Pd)

# %%
"""
# Paper permeability check
"""

# %%
Q = 0.176*1e-6/3600
A = np.pi*28e-3**2
#A = 24.63e-6
deltap = 10e6
th = 1.79e-3
P = Q*th/(deltap*A) # [m^3/msPa]
print(Q,A,deltap,th)
print(P)

# %%
import scipy.constants as sc
p1 = 10
T = 55+273.15
R = sc.R
print(p1,T,R)
Pmol = P*(p1/(T*R))
print(Pmol)

# %%
J = Pmol*deltap/th
print(J)
print(J*A)

# %%
contained = deltap/(T*R)
print(contained)
print(contained/(J*A))
print(contained/(J*A)/3600/24/1e6/365)

# %%
"""
# Arc case
"""

# %%
perm = arc(layers=3,r=0.1,extent=0.1,L=np.array([0.4,0.3,0.3]),D=np.array([1.0,0.1,0.5]),S=np.array([1.0,1.1,0.7]),\
              C0=1.0,C1=0.0,touts=[0.001,0.05,0.2,2.0],tstep=0.001,ncpu=1,N=[40,36,32])

# %%
perm.submit_job()

# %%
perm.read_field('C')

# %%
print(perm.field)

# %%
x = perm.field['C'][-1]['x']
y = perm.field['C'][-1]['y']
C = np.ravel(perm.field['C'][-1]['data'])
plt.scatter(x,y,c=C)
plt.colorbar()
plt.show()

# %%
