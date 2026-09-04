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
dir='../cont/cont_slab'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
outname='cont_chn_distr'
nhi=100

""""""
x=np.zeros((nmd,nT,nhi))
p=np.zeros((nmd,nT,nhi))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npz')['Pchn']
        P=q[q>0.0]
        hist,edge=np.histogram(np.log10(P),bins=nhi)
        x[i,j]=(edge[:-1]+edge[1:])/2
        p[i,j]=hist/np.sum(hist*(np.max(x[i,j])-np.min(x[i,j]))/nhi)

""""""
np.savez(f'{outname}.npz',x=x,p=p)