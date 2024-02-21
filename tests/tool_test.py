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
from andvaranaut import *

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
# Homogenisation Module Tests
"""

# %%
"""
## Microstructure
"""

# %%
perm = homogenisation(materials=2,D=np.array([1.0,0.1]),S=np.array([1.0,1.1]),vFrac=np.array([0.4,0.6]),\
              C0=1.0,C1=0.0,touts=[0.001,0.05,0.2,2.0,4.0],tstep=0.001,ncpu=8)

# %%
perm.reuss_mesh(Nx=1,Ny=80)

# %%
perm.submit_job()

# %%
perm.read_field('J')

# %%
perm.read_field('C')

# %%
perm.read_field('p')

# %%
print(np.mean(perm.field['J'][-1]['data'][:,1]))

# %%
perm.reuss_bound()

# %%
print(perm.D_eff,perm.S_eff,perm.P_eff)

# %%
C = perm.field['C'][-1]['data']
p = perm.field['p'][-1]['data']
Cx = perm.field['C'][-1]['x']
Cy = perm.field['C'][-1]['y']
px = perm.field['p'][-1]['x']
py = perm.field['p'][-1]['y']
Jy = perm.field['J'][-1]['data'][:,1]
mats = perm.field['J'][-1]['material']

# %%
gradC = np.zeros_like(Jy)
for i in range(len(gradC)):
    gradC[i] = Jy[i]/-perm.D[mats[i]]

# %%
print(np.mean(gradC))
print(np.mean(Jy)/(perm.C0-perm.C1))
print(np.mean(Jy)/-np.mean(gradC))
print(np.mean(gradC)/(perm.C1-perm.C0))

# %%
print(len(C),len(p),len(Jy))

# %%
print(perm.field['C'][-1]['node'])

# %%
a = np.random.randint(100,size=5)
print(a)
print(np.isin(1,a))

# %%
print([i for i in zip(C[:,0],perm.field['C'][-1]['node'])])

# %%
print(perm.field['C'][-1]['node'])

# %%
print(perm.nodesets)

# %%
g = GPMCMC(nx=2,ny=1,target=np.exp,priors=[st.norm(),st.norm()],verbose=True,kernel='RBF',noise=False)

# %%
g.set_data(np.c_[px,py],p)

# %%
g.change_conrevs([maxmin(g.x[:,i]) for i in range(2)],[meanstd(g.y[:,0])])
#g.change_conrevs([maxmin(g.x[:,0]),wgp(['maxmin','kumaraswamy'],np.ones(2),g.x[:,1])],[wgp(['stdshift','sinharcsinh','meanstd'],np.ones(3),g.y[:,0])])
#g.change_conrevs([wgp(['maxmin','kumaraswamy'],np.ones(2),g.x[:,i]) for i in range(2)],[wgp(['stdshift','sinharcsinh','meanstd'],np.ones(3),g.y[:,0])])
#g.change_conrevs([maxmin(g.x[:,i]) for i in range(2)],[wgp(['stdshift','sal','sinharcsinh','meanstd'],np.ones(7),g.y[:,0])])
g.change_model('Exponential',noise=False)

# %%
g.fit(iwgp=False,cwgp=False)

# %%
print(g.hypers)

# %%
lenC = len(C)
lenp = len(p)
perm.field['C'][-1]['S'] = np.zeros(lenC)
perm.field['C'][-1]['grad'] = np.zeros(lenC)
eps = 1e-3
predy0 = py-eps
predy1 = py+eps
gradp = (g.predict(np.c_[px,py-eps])-g.predict(np.c_[px,py+eps]))/(2*eps)
perm.field['p'][-1]['grad'] = copy.deepcopy(gradp)
deletions = []
for i in range(lenC):
    for j in range(lenp):
        if np.isclose(Cx[i],px[j]) and np.isclose(Cy[i],py[j]):
            perm.field['C'][-1]['S'][i] = C[i]/p[j]
            perm.field['C'][-1]['grad'][i] = gradp[j]*perm.field['C'][-1]['S'][i]
        if np.isclose(Cy[i],0.0) or np.isclose(Cy[i],1.0):
            deletions.append(i)
for i in range(lenC):
    for j in range(lenC):
        if i != j and np.isclose(Cx[i],Cx[j]) and np.isclose(Cy[i],Cy[j]):
            deletions.append(i)

# %%
plt.scatter(perm.field['C'][-1]['y'],perm.field['C'][-1]['data'])
plt.scatter(perm.field['C'][-1]['y'],perm.field['C'][-1]['S'])
plt.scatter(perm.field['C'][-1]['y'],perm.field['C'][-1]['grad'])
plt.show()

# %%
plt.scatter(perm.field['p'][-1]['y'],perm.field['p'][-1]['data'])
plt.scatter(perm.field['p'][-1]['y'],perm.field['p'][-1]['grad'])
plt.ylim(-2,2)
plt.show()

# %%
Cygrad = np.delete(perm.field['C'][-1]['y'],deletions,axis=0)
gradC = np.delete(perm.field['C'][-1]['grad'],deletions,axis=0)

# %%
print(gradC)

# %%
plt.scatter(perm.field['C'][-1]['y'],perm.field['C'][-1]['data'])
plt.scatter(perm.field['C'][-1]['y'],perm.field['C'][-1]['S'])
plt.scatter(Cygrad,gradC)
plt.show()

# %%
print(np.mean(gradC))

# %%
print(perm.D_eff,perm.S_eff,perm.P_eff)

# %%
print(perm.P_eff/np.mean(gradC))

# %%
0.15625/0.154

# %%
perm.cross_section_mesh(nc=15,algorithm='LS')

# %%
perm.submit_job()

# %%
perm.read_field('J')

# %%
# Volume average flux

# %%
"""
# Permeatus 2-layer planar test
"""

# %%
perm = infrastructure(materials=2,L=np.array([0.5,0.5]),D=np.array([1.0,0.1]),S=np.array([1.0,1.1]),\
              C0=1.0,C1=1e-6,touts=[0.001,0.05,0.2,2.0],tstep=0.001,ncpu=8,N=[40,36])
#perm = infrastructure(layers=2,L=np.array([0.4,0.6]),D=np.array([1.0,0.1]),S=np.array([1.0,1.1]),\
#              C0=1.0,C1=1e-6,touts=[0.001,0.05,0.2,2.0],tstep=0.001,ncpu=1,N=[24,36])
#perm = planar(layers=2,L=np.array([0.5,0.5]),P=[1.0,0.11],\
#              p0=1.0,p1=0.0,touts=[0.001,0.05,0.2,2.0],tstep=0.001,ncpu=1,N=[40,36])

# %%
print(perm.ncpu)

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
#plt.savefig('abaqus2layer.pdf',bbox_inches='tight')
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
print(p)

# %%
pdiff = np.diff(p)
print(-pdiff[0]/perm.L[0],-pdiff[1]/perm.L[1],-(pdiff[0]+pdiff[1]))

# %%
print(C)
Cdiff = np.diff(C)
print(-Cdiff[0]/perm.L[0],-Cdiff[-1]/perm.L[1],-(Cdiff[0]+Cdiff[-1])/2)

# %%
perm.get_avg_coeffs()
print(perm.Davg,perm.Savg,perm.Pavg)

# %%
print(perm.C0,perm.C1)
print(perm.p0,perm.p1)

# %%
print(1/(perm.L[0]/perm.D[0]+perm.L[1]/perm.D[1]))

# %%
print(1/(perm.L[0]/perm.S[0]+perm.L[1]/perm.S[1]))

# %%
print(1/(perm.L[0]/perm.P[0]+perm.L[1]/perm.P[1]))

# %%
print((perm.L[0]/perm.S[0]+perm.L[1]/perm.S[1])/(perm.L[0]/perm.P[0]+perm.L[1]/perm.P[1]))

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
plt.xlim(0.0,None)
plt.show()

# %%
