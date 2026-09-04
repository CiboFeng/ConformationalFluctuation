"""Import Modules"""
import argparse
import numpy as np
import os
from scipy.spatial.distance import squareform,pdist
from sklearn.decomposition import PCA
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
outname='pca'
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
pca=PCA(n_components=2)
pca.fit(D_rdc)
pc_all=pca.transform(D_rdc)
perct=pca.explained_variance_ratio_

pc=[]
for i in range(len(files)):
    pc+=[pc_all[np.sum(nfrs[:i]):np.sum(nfrs[:i+1])]]

PC1=np.zeros((len(files),nhi))
P_PC1=np.zeros((len(files),nhi))
PC1_mean=np.zeros(len(files))
PC1_std=np.zeros(len(files))
PC2=np.zeros((len(files),nhi))
P_PC2=np.zeros((len(files),nhi))
PC2_mean=np.zeros(len(files))
PC2_std=np.zeros(len(files))
for i in range(len(files)):
    dPC1=(np.max(pc[i][:,0])-np.min(pc[i][:,0]))/nhi
    hist,edge=np.histogram(pc[i][:,0],bins=nhi)
    PC1[i]=(edge[:-1]+edge[1:])/2
    P_PC1[i]=hist/np.sum(hist*dPC1)
    PC1_mean[i]=np.mean(pc[i][:,0])
    PC1_std[i]=np.std(pc[i][:,0])
    dPC2=(np.max(pc[i][:,1])-np.min(pc[i][:,1]))/nhi
    hist,edge=np.histogram(pc[i][:,1],bins=nhi)
    PC2[i]=(edge[:-1]+edge[1:])/2
    P_PC2[i]=hist/np.sum(hist*dPC2)
    PC2_mean[i]=np.mean(pc[i][:,1])
    PC2_std[i]=np.std(pc[i][:,1])

P_PC=np.zeros((len(files),nhi,nhi))
for i in range(len(files)):
    dPC1=(np.max(pc[i][:,0])-np.min(pc[i][:,0]))/nhi
    dPC2=(np.max(pc[i][:,1])-np.min(pc[i][:,1]))/nhi
    hist,xedge,yedge=np.histogram2d(pc[i][:,0],pc[i][:,1],bins=nhi)
    P_PC[i]=hist/np.sum(hist*dPC1*dPC2)

""""""
pyw=open(f'{outname}.pyw','w')
pyw.write('Percentage of the 1st, and 2nd Component:\n')
pyw.write(f'{perct[0]} {perct[1]}\n\n')
for i in range(len(files)):
    pyw.write(f'The 1st, and 2nd Component: ({labels[i]})\n')
    for j in range(nfrs[i]):
        pyw.write(f'{pc[i][j,0]} {pc[i][j,1]}\n')
    pyw.write('\n')
pyw.close()

np.savez(f'{outname}_distr.npz',PC1=PC1,P_PC1=P_PC1,PC1_mean=PC1_mean,PC1_std=PC1_std,
         PC2=PC2,P_PC2=P_PC2,PC2_mean=PC2_mean,PC2_std=PC2_std,P_PC=P_PC)