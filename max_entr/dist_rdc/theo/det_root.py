"""Import Modules"""
import argparse
import numpy as np
import os
from scipy.spatial.distance import squareform,pdist
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.figure import figaspect
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from xtc import xtc_rd

"""Set Arguments"""
xtcfile='sim_fus_0.xtc'
nat=163
ach=[0,81,162]
# ach=[0,54,108,162]

""""""
r=xtc_rd(f'{xtcfile}',f'{xtcfile[:-4]}.dat')[0]
r=r[::10]
nfr=np.shape(r)[0]
D=np.zeros((nfr,nat,nat))
for i in range(nfr):
    D[i]=squareform(pdist(r[i],'euclidean'))
    
A=np.ones((nfr,len(ach)+2,len(ach)+2))
B=np.ones((nfr,len(ach)+1,len(ach)+1))
for i in range(len(ach)):
    for j in range(len(ach)):
        A[:,i,j]=D[:,ach[i],ach[j]]**2
        B[:,i,j]=D[:,ach[i],ach[j]]**2
A[:,-1,-1]=np.zeros(nfr)
A[:,-2,-2]=np.zeros(nfr)
B[:,-1,-1]=np.zeros(nfr)

i=0
j=10
k=30
a=A[i].copy()
for l in range(len(ach)):
    a[l,len(ach)]=D[i,ach[l],j]**2
    a[len(ach),l]=D[i,ach[l],k]**2
root1=np.sqrt(-np.linalg.det(a)/np.linalg.det(B[i]))

ach_=ach+[j,k]
C=np.ones((len(ach_)+1,len(ach_)+1))
for l in range(len(ach_)):
    for m in range(len(ach_)):
        C[l,m]=D[i,ach_[l],ach_[m]]**2
C[-1,-1]=0
print('d^2_real=',C[-2,-3])
print('det(C_real)=',np.linalg.det(C))

C_=C.copy()
C_[-2,-3]=root1**2
C_[-3,-2]=root1**2
print('d^2_pred=',root1**2)
print('det(C_pred)=',np.linalg.det(C_))

x=np.linspace(min(C[-2,-3],2*root1**2-C[-2,-3])-10,max(C[-2,-3],2*root1**2-C[-2,-3])+10,10000)
y=np.zeros_like(x)
C__=C.copy()
for l in range(len(x)):
    C__[-2,-3]=x[l]
    C__[-3,-2]=x[l]
    y[l]=np.linalg.det(C__)

""""""
plt.plot(x,y,'k')
plt.plot([np.min(x),np.max(x)],[0,0],'k:')
plt.plot([C[-2,-3],C[-2,-3]],[np.min(y),np.max(y)],label=f'Real')
plt.plot([root1**2,root1**2],[np.min(y),np.max(y)],label=f'Pred')
plt.xlabel('$d_{ij}^2$')
plt.ylabel('det')
plt.legend(loc='best')
plt.savefig(f'det_root_{len(ach)}.png')
plt.show()

