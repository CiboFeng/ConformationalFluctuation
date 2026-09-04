"""Import Modules"""
import numpy as np
from scipy.optimize import curve_fit
import random
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap,ListedColormap
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from pyw import pyw

"""Set Arguments"""
dir='autocor_dist_in'
outname='proc_autocor_dist_in'
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
ac=np.zeros((nmd,nT,nfr,nat,nat))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npy')
        nfr_=np.shape(q)[0]
        ac[i,j,:nfr_]=q
ac[ac==0.0]=np.nan

"""Interpolation"""
ac_tgt=np.arange(0.1,1.0,0.1)
t_relx=np.zeros((nmd,nT,nat,nat,len(ac_tgt)))
t_relxe=np.zeros((nmd,nT,nat,nat))
for i in range(nmd):
    for j in range(nT):
        for k in range(nat):
            for l in range(k):
                for m in range(len(ac_tgt)):
                    for n in range(nfr-1):
                        if ac[i,j,n,k,l]>ac_tgt[m] and ac[i,j,n+1,k,l]<=ac_tgt[m]:
                            t_relx[i,j,k,l,m]=((ac[i,j,n+1,k,l]-ac_tgt[m])*t[n]+(ac_tgt[m]-ac[i,j,n,k,l])*t[n+1])/(ac[i,j,n+1,k,l]-ac[i,j,n,k,l])
                            t_relx[i,j,l,k,m]=t_relx[i,j,k,l,m]
                            break
                for n in range(nfr-1):
                    if ac[i,j,n,k,l]>1/np.e and ac[i,j,n+1,k,l]<=1/np.e:
                        t_relxe[i,j,k,l]=((ac[i,j,n+1,k,l]-1/np.e)*t[n]+(1/np.e-ac[i,j,n,k,l])*t[n+1])/(ac[i,j,n+1,k,l]-ac[i,j,n,k,l])
                        t_relxe[i,j,l,k]=t_relxe[i,j,k,l]
                        break

ac_mean=np.zeros((nmd,nT,nfr,nat))
t_relx_mean=np.zeros((nmd,nT,nat,len(ac_tgt)))
t_relxe_mean=np.zeros((nmd,nT,nat))
for i in range(nat):
    for j in range(nat-i):
        ac_mean[:,:,:,i]+=ac[:,:,:,j,i+j]
        t_relx_mean[:,:,i]+=t_relx[:,:,j,i+j]
        t_relxe_mean[:,:,i]+=t_relxe[:,:,j,i+j]
    ac_mean[:,:,:,i]/=nat-i
    t_relx_mean[:,:,i]/=nat-i
    t_relxe_mean[:,:,i]/=nat-i

"""Output Data"""
np.savez(f'{outname}.npz',ac_mean=ac_mean,t_relx_mean=t_relx_mean,t_relxe_mean=t_relxe_mean)
