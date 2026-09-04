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
outname='dist_rdc'
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

D_repr=np.zeros((nfr,nat,nat))
for i in range(len(ach)):
    D_repr[:,ach[i]]=D[:,ach[i]]
    D_repr[:,:,ach[i]]=D[:,:,ach[i]]
    
A=np.ones((nfr,len(ach)+2,len(ach)+2))
B=np.ones((nfr,len(ach)+1,len(ach)+1))
for i in range(len(ach)):
    for j in range(len(ach)):
        A[:,i,j]=D[:,ach[i],ach[j]]**2
        B[:,i,j]=D[:,ach[i],ach[j]]**2
A[:,-1,-1]=np.zeros(nfr)
A[:,-2,-2]=np.zeros(nfr)
B[:,-1,-1]=np.zeros(nfr)

for i in range(nfr):
    for j in range(nat):
        for k in range(j):
            if (j not in ach) and (k not in ach):
                a=A[i].copy()
                for l in range(len(ach)):
                    a[l,len(ach)]=D[i,ach[l],j]**2
                    a[len(ach),l]=D[i,ach[l],k]**2
                D_repr[i,j,k]=np.sqrt(-np.linalg.det(a)/np.linalg.det(B[i]))
                D_repr[i,k,j]=D_repr[i,j,k]

""""""
np.savez(f'{outname}_{len(ach)}.npz',D=D,D_repr=D_repr)
