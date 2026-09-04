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

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-d',type=str,help='The directory to trajectory files')
args=parser.parse_args()
xtcdir=args.d
idx='_'.join(xtcdir.split('_')[-3:])
outname=f'shp/{idx}'
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
nhi=50

""""""
seq=[res.index(seq[i]) for i in range(nat)]
m=np.array([m[seq[i]] for i in range(nat)])

r=np.zeros((0,nat,3))
xtcfiles=[f for f in os.listdir(xtcdir) if f.endswith('.xtc')]
for i in range(len(xtcfiles)):
    try:
        r=np.concatenate((r,xtc_rd(f'{xtcdir}/{xtcfiles[i]}',f'{xtcdir}/{xtcfiles[i][:-4]}.dat')[0]),axis=0)
    except Exception as err:
        print(err)

nfr=np.shape(r)[0]
shp=np.zeros((6,nfr))
r-=np.sum(m[np.newaxis,:,np.newaxis]*r,axis=1,keepdims=True)/np.sum(m)
IT=np.sum(np.einsum('ijk,ijl->ijkl',np.sqrt(m[np.newaxis,:,np.newaxis])*r,np.sqrt(m[np.newaxis,:,np.newaxis])*r),axis=1)/np.sum(m)
tr=np.trace(IT,axis1=1,axis2=2)
lmd=np.array([np.sort(np.linalg.eigvals(IT[i]))[::-1] for i in range(nfr)])
shp[:3,:]=(np.sqrt(lmd)).T
shp[3,:]=np.sqrt(tr)
shp[4,:]=9/2*np.var(lmd,axis=1)/tr**2
shp[5,:]=27*np.product(lmd-np.mean(lmd,axis=1,keepdims=True),axis=1)/tr**3

shp_mean=np.mean(shp,axis=1)
shp_std=np.std(shp,axis=1)

Shp=np.zeros((6,nhi))
P1=np.zeros((6,nhi))
for i in range(6):
    dShp=(np.max(shp[i])-np.min(shp[i]))/nhi
    hist,edge=np.histogram(shp[i],bins=nhi)
    Shp[i]=(edge[:-1]+edge[1:])/2
    P1[i]=hist/np.sum(hist*dShp)

P2=np.zeros((6,6,nhi,nhi))
for i in range(6):
    for j in range(6):
        if i>j:
            dShp1=(np.max(shp[i])-np.min(shp[i]))/nhi
            dShp2=(np.max(shp[j])-np.min(shp[j]))/nhi
            hist,xedge,yedge=np.histogram2d(shp[i],shp[j],bins=nhi)
            P2[i,j]=hist/np.sum(hist*dShp1*dShp2)
            P2[j,i]=P2[i,j].T

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.savez(f'{outname}.npz',shp=shp,shp_mean=shp_mean,shp_std=shp_std,Shp=Shp,P1=P1,P2=P2)