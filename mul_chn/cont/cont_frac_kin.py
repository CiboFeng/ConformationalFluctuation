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
dir0='cont_slab'
dir='cont_kin'
mdls=['1.37sgm_0.6sgm','1.5sgm_0.82sgm','1.55sgm_1.0sgm']
nmd=len(mdls)
Mdls=['Flex','Mid','Fix']
outname='cont_frac_kin'
nfr=10000
nat=163
nrep=20
dt=10e-6*10000 ### ns

""""""
# Pintra0=np.zeros((nmd,nat,nat))
# Pinter0=np.zeros((nmd,nat,nat))
# for i in range(nmd):
#     q=np.load(f'{dir0}/{mdls[i]}_300.0.npz')
#     Pintra0[i]=q['Pintra']
#     Pinter0[i]=q['Pinter']
Pintra=np.zeros((nmd,nfr,nat,nat))
Pinter=np.zeros((nmd,nfr,nat,nat))
Qintra=np.zeros((nmd,nfr))
Qinter=np.zeros((nmd,nfr))
for i in range(nmd):
    for j in range(nrep):
        q=np.load(f'{dir}/{mdls[i]}_{j}.npz')
        Pintra[i]+=q['Pintra']
        Pinter[i]+=q['Pinter']
    Pintra[i]/=nrep
    Pinter[i]/=nrep
    for j in range(nfr):
        # Qintra[i,j]=r2_score(Pintra[i,j][np.triu_indices(nat,k=4)],Pintra0[i][np.triu_indices(nat,k=4)])
        # Qinter[i,j]=r2_score(Pinter[i,j].reshape(-1),Pinter0[i].reshape(-1))
        Qintra[i,j]=np.sum(Pintra[i,j][np.triu_indices(nat,k=4)])
        Qinter[i,j]=np.sum(Pinter[i,j])

""""""
pyw=open(f'{outname}.pyw','w')
for i in range(nmd):
    pyw.write(f'Time (ns), Intra-chain, and Inter-chain Contact Fraction: ({Mdls[i]})\n')
    for j in range(nfr):
        pyw.write(f'{j*dt} {Qintra[i,j]} {Qinter[i,j]}\n')
    pyw.write('\n')
pyw.close()
