"""Import Modules"""
import argparse
import numpy as np
import os
from scipy.spatial.distance import squareform,pdist
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.figure import figaspect
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from pyw import pyw

"""Set Arguments"""
dir1='cont_slab'
dir2='cont_kin'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
outname='cont_kin_show'
nfr=10000
nch1=100
nch2=50
nat=163
nrep=20
dt=10e-6*10000 ### ns

""""""
Pintra1=np.zeros((nmd,nT,nat,nat))
Pinter1=np.zeros((nmd,nT,nat,nat))
for i in range(nmd):
    for j in range(nT):
        if j==0:
            q=np.load(f'{dir1}/{mdls[i]}_{Ts[j]}.npz')
            Pintra1[i,j]=q['Pintra']
            Pinter1[i,j]=q['Pinter']
Pintra2=np.zeros((nmd,nfr,nat,nat))
Pinter2=np.zeros((nmd,nfr,nat,nat))
for i in range(nmd):
    for j in range(nrep):
        q=np.load(f'{dir2}/{mdls[i]}_{j}.npz')
        Pintra2[i]+=q['Pintra']
        Pinter2[i]+=q['Pinter']
    Pintra2[i]/=nrep
    Pinter2[i]/=nrep

idx=[0,10,100,1000,9999]
Pintra2=Pintra2[:,idx]
Pinter2=Pinter2[:,idx]

""""""
np.savez(f'{outname}.npz',Pintra1=Pintra1,Pinter1=Pinter1,Pintra2=Pintra2,Pinter2=Pinter2)