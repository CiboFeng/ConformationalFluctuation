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
outname=f'dist/{mdl}_{T}'
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]
nhiz=250
nhid=50
dlim=[0.0,80.0]

""""""
m=np.array([m[seq[i]] for i in range(nat)])
r,tp,b=xtc_rd(xtcfile,f'{xtcfile[:-18]}.dat')
nfr,N,_=np.shape(r)
nch=round(N/nat)
r=r.reshape(nfr,nch,nat,3)
D=np.zeros((nfr,nch,nat,nat))
for i in range(nfr):
    for j in range(nch):
        D[i,j]=squareform(pdist(r[i,j],'euclidean'))
zc=np.sum(m[np.newaxis,np.newaxis]*r[:,:,:,2],axis=-1)/np.sum(m)

binz=np.linspace((b[0,-1,0]-b[0,-1,1])/2,(b[0,-1,1]-b[0,-1,0])/2,nhiz+1)
bind=np.linspace(dlim[0],dlim[1],nhid+1)
z=(binz[:-1]+binz[1:])/2
d=(bind[:-1]+bind[1:])/2

d2_mean=np.zeros((nat,nat,nhiz))
d2_std=np.zeros((nat,nat,nhiz))
for i in range(nhiz):
    D_=D[(zc>binz[i]) & (zc<=binz[i+1])] ### Can't be tht[binz[i]<zc<=binz[i+1]]
    if np.shape(D_)[0]>0:
        d2_mean[:,:,i]=np.mean(D_,axis=0)
        d2_std[:,:,i]=np.std(D_,axis=0)
p2=np.zeros((nat,nat,nhiz,nhid))
p1=np.zeros((nat,nhiz,nhid))
d1_mean=np.zeros((nat,nhiz))
d1_std=np.zeros((nat,nhiz))
for i in range(nat):
    for j in range(i):
        hist,xedge,yedge=np.histogram2d(zc.reshape(-1),D[:,:,i,j].reshape(-1),bins=[binz,bind])
        p2[i,j]=hist/np.sum(hist*(np.max(z)-np.min(z))/nhiz*(np.max(d)-np.min(d))/nhid)
        p2[j,i]=p2[i,j]
    D_tmp=np.zeros((nfr,nch,nat-i))
    for j in range(nat-i):
        D_tmp[:,:,j]=D[:,:,j,i+j]
    Zc=np.tile(zc[:,:,np.newaxis],(1,1,nat-i))
    hist,xedge,yedge=np.histogram2d(Zc.reshape(-1),D_tmp.reshape(-1),bins=[binz,bind])
    p1[i]=hist/np.sum(hist*(np.max(z)-np.min(z))/nhiz*(np.max(d)-np.min(d))/nhid)
    for j in range(nhiz):
        D_=D_tmp[(zc>binz[j]) & (zc<=binz[j+1])] ### Can't be tht[binz[i]<zc<=binz[i+1]]
        if np.shape(D_)[0]>0:
            d1_mean[i,j]=np.mean(D_)
            d1_std[i,j]=np.std(D_)

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.savez(f'{outname}.npz',z=z,d=d,p2=p2,d2_mean=d2_mean,d2_std=d2_std,p1=p1,d1_mean=d1_mean,d1_std=d1_std)