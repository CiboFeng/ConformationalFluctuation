"""Import Modules"""
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

"""Set Arguments"""
pywfile='pca.pyw'
npzfiles=['../../all_atom/shp/shp.npz',
          '../../hps/shp/shp.npz',
         '../shp/shp/1.55sgm_1.0sgm_itr20.npz',
         '../shp/shp/1.5sgm_0.82sgm_itr20.npz',
         '../shp/shp/1.37sgm_0.6sgm_itr20.npz']
outname='pca_shp'
labels=['AA','HPS','Qch','Mid','Flx']
nmd=len(labels)-2
nhi=50

"""Read Data and Calculate"""
pc=[]*(nmd+2)
for i in range(nmd+2):
    pc+=[pyw(pywfile,f'({labels[i]})')]

shp=[]
for i in range(nmd+2):
    q=np.load(npzfiles[i])
    shp+=[q['shp']]

q=[]
for i in range(nmd+2):
    q+=[np.concatenate((pc[i],shp[i]),axis=0)]
nq=np.shape(q[0])[0]

cor=np.zeros((nmd+2,nq,nq))
q_mean=np.zeros((nmd+2,nq))
q_std=np.zeros((nmd+2,nq))
Q=np.zeros((nmd+2,nq,nhi))
P1=np.zeros((nmd+2,nq,nhi))
P2=np.zeros((nmd+2,nq,nq,nhi,nhi))
for i in range(nmd+2):
    cor[i]=np.corrcoef(q[i])

    q_mean[i]=np.mean(q[i],axis=1)
    q_std[i]=np.std(q[i],axis=1)

    for j in range(nq):
        dQ=(np.max(q[i][j])-np.min(q[i][j]))/nhi
        hist,edge=np.histogram(q[i][j],bins=nhi)
        Q[i,j]=(edge[:-1]+edge[1:])/2
        P1[i,j]=hist/np.sum(hist*dQ)
        for k in range(nq):
            if j>k:
                dQ1=(np.max(q[i][j])-np.min(q[i][j]))/nhi
                dQ2=(np.max(q[i][k])-np.min(q[i][k]))/nhi
                hist,xedge,yedge=np.histogram2d(q[i][j],q[i][k],bins=nhi)
                P2[i,j,k]=hist/np.sum(hist*dQ1*dQ2)
                P2[i,k,j]=P2[i,j,k].T

"""Output Data"""
np.savez(f'{outname}.npz',cor=cor,q_mean=q_mean,q_std=q_std,Q=Q,P1=P1,P2=P2)