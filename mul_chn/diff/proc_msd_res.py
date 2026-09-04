"""Import Modules"""
import numpy as np
from scipy.optimize import curve_fit
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from pyw import pyw

"""Set Arguments"""
dir='msd_res'
outname='proc_msd_res'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nfr=10000
nch=100
nat=163
dt=10e-6*10000 ### ns

"""Read Data"""
t=np.arange(nfr)*dt
d=np.zeros((nmd,nT,nfr,nch,nat,3))
for i in range(nmd):
    for j in range(nT):
        d[i,j]=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npy')

dm=np.mean(d,axis=-3)
dmm=np.mean(d,axis=(-3,-2))

"""Plot"""
np.savez(f'{outname}.npz',dm=dm,dmm=dmm)