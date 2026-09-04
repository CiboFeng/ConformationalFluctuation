"""Import Modules"""
import argparse
import numpy as np
import os
from scipy.spatial.distance import squareform,pdist
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.figure import figaspect
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from pyw import pyw
from xtc import xtc_rd

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-f',type=str,help='The trajectory file')
args=parser.parse_args()
xtcfile=args.f
mdl='_'.join(xtcfile.split('/')[-2].split('_')[-3:-1])
T=xtcfile.split('/')[-2].split('_')[-1]
pywfile=f'../dst/dst/{mdl}_{T}.pyw'
outname=f'res_distr/{mdl}_{T}'
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]
ntp=len(res)
nhi=250

"""Read Data"""
m=np.array([m[seq[i]] for i in range(nat)])
r,tp,b=xtc_rd(xtcfile,f'{xtcfile[:-9]}.dat')
nfr,N,_=np.shape(r)
nch=round(N/nat)
r=r.reshape(nfr,nch,nat,3)
z=r[:,:,:,-1]

z1,z2,d=pyw(pywfile,'Position')[:3,0]

""""""
zabs=np.abs(z-(z1+z2)/2)
zabs_mean=np.mean(zabs,axis=(0,1))
zabs_std=np.std(zabs,axis=(0,1))

bin=np.linspace((b[0,-1,0]-b[0,-1,1])/2,(b[0,-1,1]-b[0,-1,0])/2,nhi+1)
idx=np.tile(np.arange(nat),(nfr,nch,1))
hist,xedge,yedge=np.histogram2d(z.reshape(-1),idx.reshape(-1),bins=[bin,nat])
Z=(xedge[:-1]+xedge[1:])/2
pZ=hist/np.sum(hist*(np.max(Z)-np.min(Z))/nhi)

""""""
z_tp=[np.zeros((nfr,nch,0)) for _ in range(ntp)] ### [[]]*ntp will lead to the 20 sub-lists change synchronously.
for i in range(nat):
    z_tp[seq[i]]=np.concatenate((z_tp[seq[i]],z[:,:,i:i+1]),axis=-1)

zabs_tp_mean=np.zeros(ntp)
zabs_tp_std=np.zeros(ntp)
for i in range(ntp):
    zabs=np.abs(z_tp[i]-(z1+z2)/2)
    if np.shape(zabs)[-1]>0:
        zabs_tp_mean[i]=np.mean(zabs)
        zabs_tp_std[i]=np.std(zabs)

idx=np.tile(seq,(nfr,nch,1))
bin=np.linspace((b[0,-1,0]-b[0,-1,1])/2,(b[0,-1,1]-b[0,-1,0])/2,nhi+1)
bintp=np.linspace(-0.5,ntp-0.5,ntp+1)
hist,xedge,yedge=np.histogram2d(z.reshape(-1),idx.reshape(-1),bins=[bin,bintp])
Z_tp=(xedge[:-1]+xedge[1:])/2
pZ_tp=hist/np.sum(hist*(np.max(Z_tp)-np.min(Z_tp))/nhi)

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.savez(f'{outname}.npz',Z=Z,pZ=pZ,zabs_mean=zabs_mean,zabs_std=zabs_std,
         Z_tp=Z_tp,pZ_tp=pZ_tp,zabs_tp_mean=zabs_tp_mean,zabs_tp_std=zabs_tp_std)