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
outname=f'mag/{mdl}_{T}'
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]
nhi_z=250
nhi_mag=50
maglim=[[-1,1],[-1,1],[-1,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1]]
nq=len(maglim)

""""""
m=np.array([m[seq[i]] for i in range(nat)])
r,tp,b=xtc_rd(xtcfile,f'{xtcfile[:-18]}.dat')
nfr,N,_=np.shape(r)
nch=round(N/nat)
r=r.reshape(nfr,nch,nat,3)

z=r[:,:,:,-1]
zc=np.sum(m[np.newaxis,np.newaxis]*z,axis=-1)/np.sum(m)
r-=np.sum(m[np.newaxis,np.newaxis,:,np.newaxis]*r,axis=-2,keepdims=True)/np.sum(m)

rN=np.mean(r[:,:,:50],axis=-2)
rC=np.mean(r[:,:,-50:],axis=-2)
drNC=rC-rN
ev=np.zeros((nfr,nch,3))
for i in range(nch):
    IT=np.sum(np.einsum('ijk,ijl->ijkl',np.sqrt(m[np.newaxis,:,np.newaxis])*r[:,i],np.sqrt(m[np.newaxis,:,np.newaxis])*r[:,i]),axis=1)/np.sum(m)
    ev_tmp=[np.linalg.eigh(IT[j])[1][:,-1] for j in range(nfr)] ### eigh function is more suitable for hermitian matrix. The z component of the last column is cos_theta.
    ev[:,i]=ev_tmp*np.sign(np.sum(ev_tmp*drNC[:,i],axis=-1,keepdims=True))

bin_z=np.linspace((b[0,-1,0]-b[0,-1,1])/2,(b[0,-1,1]-b[0,-1,0])/2,nhi_z+1)
# def assort(x):
#     y=x-x
#     for i in range(1,nhi_z):
#         y+=np.heaviside(x-bin_z[i],0.0)
#     return y.astype(int)
mag=np.zeros((nfr,nhi_z,nq))
for i in range(nfr):
    # idx=[[] for _ in range(nhi_z)]
    # idx_bin=assort(zc[i])
    # for j in range(nch):
    #     idx[idx_bin[j]]+=[j]
    # for j in range(nhi_z):
    #     mag[i,j]=np.linalg.norm(np.sum(ev[i,idx[j]],axis=0))
    for j in range(nhi_z):
        idx=(zc[i]>bin_z[j]) & (zc[i]<=bin_z[j+1])
        ev_tmp=ev[i,idx]
        ev_mean=np.mean(ev_tmp,axis=0)
        mag[i,j,:3]=ev_mean
        mag[i,j,3:6]=np.abs(ev_mean)
        mag[i,j,6:9]=np.mean(np.abs(ev_tmp),axis=0)
        mag[i,j,9]=np.linalg.norm(ev_mean)

bin_mag=np.zeros((nhi_mag+1,nq))
for i in range(nq):
    bin_mag[:,i]=np.linspace(maglim[i][0],maglim[i][1],nhi_mag+1)

mag_mean=np.nanmean(mag,axis=0)
mag_std=np.nanstd(mag,axis=0)

z=(bin_z[:-1]+bin_z[1:])/2
Mag=(bin_mag[:-1]+bin_mag[1:])/2
p=np.zeros((nhi_z,nhi_mag,nq))
for i in range(nhi_z):
    for j in range(nq):
        hist,edge=np.histogram(mag[:,i,j],bins=bin_mag[:,j])
        p[i,:,j]=hist/np.sum(hist*(np.max(Mag[:,j])-np.min(Mag[:,j]))/nhi_mag)

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.savez(f'{outname}.npz',z=z,Mag=Mag,p=p,mag_mean=mag_mean,mag_std=mag_std)