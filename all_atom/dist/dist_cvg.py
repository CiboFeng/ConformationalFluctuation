"""Import Modules"""
import argparse
import os
import numpy as np
from scipy.spatial.distance import squareform,pdist
import sys
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
from xtc import xtc_rd

""""""
parser=argparse.ArgumentParser()
parser.add_argument('-d',type=str,help='The directory to trajectory files')
args=parser.parse_args()
xtcdir=args.d
outname='dist_cvg'
nat=163
nhist=50

""""""
r=np.zeros((0,nat,3))
xtcfiles=[f for f in os.listdir(xtcdir) if f.endswith('.xtc')]
for i in range(len(xtcfiles)):
    try:
        r=np.concatenate((r,xtc_rd(f'{xtcdir}/{xtcfiles[i]}',f'{xtcdir}/{xtcfiles[i][:-4]}.pdb')[0]),axis=0)
    except Exception as err:
        print(err)

nfr,nat=np.shape(r)[:2]
D1=np.zeros((int(nfr/2),nat,nat))
D2=np.zeros((int(nfr/2),nat,nat))
for i in range(int(nfr/2)):
    D1[i]=squareform(pdist(r[i],'euclidean'))
    D2[i]=squareform(pdist(r[-i-1],'euclidean'))

""""""
x1=np.zeros((nat,nat,nhist))
px1=np.zeros((nat,nat,nhist))
x2=np.zeros((nat,nat,nhist))
px2=np.zeros((nat,nat,nhist))
for i in range(nat):
    for j in range(i):
        hist,edge=np.histogram(D1[:,i,j],bins=nhist)
        x1[i,j]=(edge[:-1]+edge[1:])/2
        x1[j,i]=x1[i,j]
        px1[i,j]=hist/np.sum(hist*(np.max(x1[i,j])-np.min(x1[i,j]))/nhist)
        px1[j,i]=px1[i,j]
        hist,edge=np.histogram(D2[:,i,j],bins=nhist)
        x2[i,j]=(edge[:-1]+edge[1:])/2
        x2[j,i]=x2[i,j]
        px2[i,j]=hist/np.sum(hist*(np.max(x2[i,j])-np.min(x2[i,j]))/nhist)
        px2[j,i]=px2[i,j]

""""""
np.savez(f'{outname}.npz',x1=x1,px1=px1,x2=x2,px2=px2)
