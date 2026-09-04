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

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-f',type=str,help='The trajectory file')
args=parser.parse_args()
xtcfile=args.f
mdl='_'.join(xtcfile.split('/')[-2].split('_')[-3:-1])
T=xtcfile.split('/')[-2].split('_')[-1]
outname=f'rg_slab/{mdl}_{T}'
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]
nhiz=250
nhirg=50
rglim=[10.0,60.0]

""""""
m=np.array([m[seq[i]] for i in range(nat)])
r,tp,b=xtc_rd(xtcfile,f'{xtcfile[:-18]}.dat')
nfr,N,_=np.shape(r)
nch=round(N/nat)
r=r.reshape(nfr,nch,nat,3)
rc=np.sum(m[np.newaxis,np.newaxis,:,np.newaxis]*r,axis=2,keepdims=True)/np.sum(m)
rg=np.sum(m[np.newaxis,np.newaxis]*np.linalg.norm(r-rc,axis=-1),axis=-1)/np.sum(m)

z=r[:,:,:,-1]
zc=np.sum(m[np.newaxis,np.newaxis]*z,axis=-1)/np.sum(m)

binz=np.linspace((b[0,-1,0]-b[0,-1,1])/2,(b[0,-1,1]-b[0,-1,0])/2,nhiz+1)
binrg=np.linspace(rglim[0],rglim[1],nhirg+1)

rg_mean=np.zeros(nhiz)
rg_std=np.zeros(nhiz)
for i in range(nhiz):
    rg_=rg[(zc>binz[i]) & (zc<=binz[i+1])] ### Can't be tht[binz[i]<zc<=binz[i+1]]
    if len(rg_)>0:
        rg_mean[i]=np.mean(rg_)
        rg_std[i]=np.std(rg_)

hist,xedge,yedge=np.histogram2d(zc.reshape(-1),rg.reshape(-1),bins=[binz,binrg])
z=(xedge[:-1]+xedge[1:])/2
rg=(yedge[:-1]+yedge[1:])/2
p=hist/np.sum(hist*(np.max(z)-np.min(z))/nhiz*(np.max(rg)-np.min(rg))/nhirg)

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.savez(f'{outname}.npz',z=z,rg=rg,p=p,rg_mean=rg_mean,rg_std=rg_std)