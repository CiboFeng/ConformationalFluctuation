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
files=['../../all_atom/dist_rdc/dist_full.npz',
       '../../hps/dist_rdc/dist_full.npz',
       'dist_rdc/1.55sgm_1.0sgm_itr20.npy',
       'dist_rdc/1.5sgm_0.82sgm_itr20.npy',
       'dist_rdc/1.37sgm_0.6sgm_itr20.npy']
labels=['AA','HPS','Qch','Mid','Flx']
outname='pca_full'
nat=163

""""""
D=np.zeros((0,nat,nat))
nfrs=np.zeros(len(files),dtype=int)
for i in range(len(files)):
    D_tmp=np.load(files[i])['D']
    D=np.concatenate((D,D_tmp),axis=0)
    nfrs[i]=np.shape(D_tmp)[0]

D=D.reshape(np.sum(nfrs),-1)
pca=PCA(n_components=2)
pca.fit(D)
pc_all=pca.transform(D)
perct=pca.explained_variance_ratio_

pc=[]
for i in range(len(files)):
    pc+=[pc_all[np.sum(nfrs[:i]):np.sum(nfrs[:i+1])]]

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