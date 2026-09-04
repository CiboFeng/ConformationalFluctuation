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
from pyw import pyw

"""Set Arguments"""
filep='prs/rerun_prs'
outname=f'surf_tens'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nhi=20
nset=100
lz=500.0

""""""
F=np.zeros((nmd,nT,nhi))
p_F=np.zeros((nmd,nT,nhi))
F_mean=np.zeros((nmd,nT))
F_std=np.zeros((nmd,nT))
for i in range(nmd):
    for j in range(nT):
        f=[]
        log=open(f'{filep}_{mdls[i]}_{Ts[j]}.log')
        lslog=log.readlines()
        for k in range(len(lslog)):
            if lslog[k][:16]=='Step Pxx Pyy Pzz':
                for l in range(k+1,len(lslog)):
                    ls_log=lslog[l].strip('\n').split()
                    try:
                        float(ls_log[0])
                    except Exception:
                        continue
                    Pxx=float(ls_log[1])
                    Pyy=float(ls_log[2])
                    Pzz=float(ls_log[3])
                    f+=[lz/2*(Pzz-(Pxx+Pyy)/2)]
        # lenset=int(np.floor(len(f)/nset))
        # f=[np.mean(f[k*lenset:(k+1)*lenset]) for k in range(nset)]
        F_mean[i,j]=np.mean(f)
        F_std[i,j]=np.std(f)
        hist,edge=np.histogram(f,bins=nhi)
        F[i,j]=(edge[:-1]+edge[1:])/2
        p_F[i,j]=hist/np.sum(hist*(np.max(F[i,j])-np.min(F[i,j]))/nhi)

""""""
np.savez(f'{outname}.npz',F=F,p_F=p_F,F_mean=F_mean,F_std=F_std)
                    