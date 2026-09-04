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
files0=['cont_slab/1.37sgm_0.6sgm_300.0.pyw',
        'cont_slab/1.5sgm_0.82sgm_300.0.pyw',
        'cont_slab/1.55sgm_1.0sgm_300.0.pyw']
mdls=['cont_kin/1.37sgm_0.6sgm',
      'cont_kin/1.5sgm_0.82sgm',
      'cont_kin/1.55sgm_1.0sgm']
idxs=['Flex','Mid','Fix']
outname='cont_frac_kin'
nat=163
nrep=20
dt=10e-6*10000 ### ns

""""""
Pintra0=np.zeros((len(mdls),nat,nat))
Pinter0=np.zeros((len(mdls),nat,nat))
for i in range(len(mdls)):
    q=pyw(files0[i],'Intra-chain,').reshape(2,nat,nat)
    Pintra0[i]=q[0]
    Pinter0[i]=q[1]
###
nfrs=np.zeros((len(mdls),nrep))
for i in range(len(mdls)):
    for j in range(nrep):
        q=np.load(f'{mdls[i]}_{j}.npz')
        nfrs[i,j]=np.shape(q['Pintra'])[0]
nfrmax=int(np.max(nfrs))
###
Pintra=np.zeros((len(mdls),nfrmax,nat,nat))
Pinter=np.zeros((len(mdls),nfrmax,nat,nat))
Qintra=np.zeros((len(mdls),nfrmax))
Qinter=np.zeros((len(mdls),nfrmax))
num=np.zeros((len(mdls),nfrmax),dtype=int)
for i in range(len(mdls)):
    for j in range(nrep):
        q=np.load(f'{mdls[i]}_{j}.npz')
        nfr=np.shape(q['Pintra'])[0]
        num[i,:nfr]+=1
        Pintra[i,:nfr]+=q['Pintra']
        Pinter[i,:nfr]+=q['Pinter']
    for j in range(nfrmax):
        try:
            Pintra[i,j]/=num[i,j]
            Pinter[i,j]/=num[i,j]
            Qintra[i,j]=r2_score(Pintra[i,j][np.triu_indices(nat,k=4)],Pintra0[i][np.triu_indices(nat,k=4)])
            Qinter[i,j]=r2_score(Pinter[i,j].reshape(-1),Pinter0[i].reshape(-1))
        except:
            pass

""""""
pyw=open(f'{outname}.pyw','w')
for i in range(len(mdls)):
    pyw.write(f'Time (ns), Intra-chain, and Inter-chain Contact Fraction: ({idxs[i]})\n')
    for j in range(nfrmax):
        pyw.write(f'{j*dt} {Qintra[i,j]} {Qinter[i,j]}\n')
    pyw.write('\n')
pyw.close()
