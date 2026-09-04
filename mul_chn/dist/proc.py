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
dir='dist'
outname='proc'
mdls=['1.55sgm_1.0sgm','1.5sgm_0.82sgm','1.37sgm_0.6sgm']
nmd=len(mdls)
Ts=[300.0,350.0,400.0,450.0]
nT=len(Ts)
nat=163
nhiz=250
nhid=50

"""Read Data"""
z=np.load(f'{dir}/{mdls[0]}_{Ts[0]}.npz')['z']
p1=np.zeros((nmd,nT,nat,nhiz,nhid))
d1_mean=np.zeros((nmd,nT,nat,nhiz))
for i in range(nmd):
    for j in range(nT):
        q=np.load(f'{dir}/{mdls[i]}_{Ts[j]}.npz')
        p1[i,j]=q['p1']
        d1_mean[i,j]=q['d1_mean']
pz=np.sum(p1,axis=-1)
pz/=np.sum(pz*(np.max(z)-np.min(z))/nhiz,axis=-1,keepdims=True)

"""Output Data"""
np.savez(f'{outname}.npz',z=z,pz=pz,d1_mean=d1_mean)
