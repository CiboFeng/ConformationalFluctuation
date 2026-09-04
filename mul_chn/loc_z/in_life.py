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

""""""
pywdir='dst'
mdls=['1.37sgm_0.6sgm','1.5sgm_0.82sgm','1.55sgm_1.0sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0]
nT=len(Ts)
nch=100

""""""
for i in range(nmd):
    for j in range(nT):
        inout=pyw(f'{pywdir}/{mdls[i]}_{Ts[j]}.pyw','In (1)')[2]
        nfr=round(len(inout)/nch)
        inout=inout.reshape(nfr,nch)
        sl=[]
        lensl=np.zeros(nch,dtype=int)
        for k in range(nch):
            diff=inout[1:,k]-inout[:-1,k]
            if np.sum(diff)==0.0:
                diff=np.concatenate((np.array([1.0]),diff))
                diff=np.concatenate((diff,np.array([-1.0])))
            else:
                for l in range(nfr-1):
                    if diff[l]!=0.0:
                        if diff[l]==-1.0:
                            diff=np.concatenate((np.array([1.0]),diff))
                        else:
                            diff=np.concatenate((np.array([0.0]),diff))
                        break
                for l in range(nfr-1,-1,-1):
                    if diff[l]!=0.0:
                        if diff[l]==1.0:
                            diff=np.concatenate((diff,np.array([-1.0])))
                        else:
                            diff=np.concatenate((diff,np.array([0.0])))
                        break
            sl.append([])
            for l in range(nfr+1):
                if diff[l]==1.0:
                    sl[-1].append([l])
                if diff[l]==-1.0:
                    sl[-1][-1]+=[l]
            lensl[k]=len(sl[k])
        print(i,j,np.min(lensl),np.mean(lensl),np.max(lensl))