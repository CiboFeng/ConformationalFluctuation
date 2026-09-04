"""Import Modules"""
import argparse
import numpy as np
import umap
from scipy.spatial.distance import squareform,pdist
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from xtc import xtc_rd

"""Set Arguments"""
files=['../../all_atom/dist_rdc/dist_rdc.npy',
       '../../hps/dist_rdc/dist_rdc.npy',
       'dist_rdc/1.55sgm_1.0sgm_itr20.npy',
       'dist_rdc/1.5sgm_0.82sgm_itr20.npy',
       'dist_rdc/1.37sgm_0.6sgm_itr20.npy']
labels=['AA','HPS','Qch','Mid','Flx']
outname='umapa'
nat=163
ach=[0,54,108,162]
nhi=50

""""""
D_rdc=np.zeros((0,len(ach),nat))
nfrs=np.zeros(len(files),dtype=int)
for i in range(len(files)):
    D_tmp=np.load(files[i])
    D_rdc=np.concatenate((D_rdc,D_tmp),axis=0)
    nfrs[i]=np.shape(D_tmp)[0]

D_rdc=D_rdc.reshape(np.sum(nfrs),-1)
reducer=umap.UMAP(n_neighbors=15,n_components=2,min_dist=1.0,metric='euclidean')
umap_all=reducer.fit_transform(D_rdc)

umap=[]
for i in range(len(files)):
    umap+=[umap_all[np.sum(nfrs[:i]):np.sum(nfrs[:i+1])]]

UMAP1=np.zeros((len(files),nhi))
P_UMAP1=np.zeros((len(files),nhi))
UMAP1_mean=np.zeros(len(files))
UMAP1_std=np.zeros(len(files))
UMAP2=np.zeros((len(files),nhi))
P_UMAP2=np.zeros((len(files),nhi))
UMAP2_mean=np.zeros(len(files))
UMAP2_std=np.zeros(len(files))
for i in range(len(files)):
    dUMAP1=(np.max(umap[i][:,0])-np.min(umap[i][:,0]))/nhi
    hist,edge=np.histogram(umap[i][:,0],bins=nhi)
    UMAP1[i]=(edge[:-1]+edge[1:])/2
    P_UMAP1[i]=hist/np.sum(hist*dUMAP1)
    UMAP1_mean[i]=np.mean(umap[i][:,0])
    UMAP1_std[i]=np.std(umap[i][:,0])
    dUMAP2=(np.max(umap[i][:,1])-np.min(umap[i][:,1]))/nhi
    hist,edge=np.histogram(umap[i][:,1],bins=nhi)
    UMAP2[i]=(edge[:-1]+edge[1:])/2
    P_UMAP2[i]=hist/np.sum(hist*dUMAP2)
    UMAP2_mean[i]=np.mean(umap[i][:,1])
    UMAP2_std[i]=np.std(umap[i][:,1])

P_UMAP=np.zeros((len(files),nhi,nhi))
for i in range(len(files)):
    dUMAP1=(np.max(umap[i][:,0])-np.min(umap[i][:,0]))/nhi
    dUMAP2=(np.max(umap[i][:,1])-np.min(umap[i][:,1]))/nhi
    hist,xedge,yedge=np.histogram2d(umap[i][:,0],umap[i][:,1],bins=nhi)
    P_UMAP[i]=hist/np.sum(hist*dUMAP1*dUMAP2)

""""""
pyw=open(f'{outname}.pyw','w')
for i in range(len(files)):
    pyw.write(f'The 1st, and 2nd Component: ({labels[i]})\n')
    for j in range(nfrs[i]):
        pyw.write(f'{umap[i][j,0]} {umap[i][j,1]}\n')
    pyw.write('\n')
pyw.close()

np.savez(f'{outname}_distr.npz',UMAP1=UMAP1,P_UMAP1=P_UMAP1,UMAP1_mean=UMAP1_mean,UMAP1_std=UMAP1_std,
         UMAP2=UMAP2,P_UMAP2=P_UMAP2,UMAP2_mean=UMAP2_mean,UMAP2_std=UMAP2_std,P_UMAP=P_UMAP)