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
mdl='_'.join(xtcfile.split('/')[-2].split('_')[-2:])
nrep=20
outname=f'mag_r/{mdl}'
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]
nhi_r=250
nhi_mag=50

""""""
m=np.array([m[seq[i]] for i in range(nat)])
r=[]
b=[]
for i in range(nrep):
    xtcfile_=xtcfile.replace('__',f'_{i}_')
    r_,tp,b_=xtc_rd(xtcfile_,f'{xtcfile_[:-9]}.dat')
    r+=[r_]
    b+=[b_]
r=np.concatenate(r,axis=0)
b=np.concatenate(b,axis=0)
nfr,N,_=np.shape(r)
nch=round(N/nat)
r=r.reshape(nfr,nch,nat,3)

rc=np.sum(m[np.newaxis,np.newaxis,:,np.newaxis]*r,axis=-2)/np.sum(m)
r-=rc[:,:,np.newaxis,:]

rN=np.mean(r[:,:,:50],axis=-2)
rC=np.mean(r[:,:,-50:],axis=-2)
drNC=rC-rN
ev=np.zeros((nfr,nch,3))
evr=np.zeros((nfr,nch))
for i in range(nch):
    IT=np.sum(np.einsum('ijk,ijl->ijkl',np.sqrt(m[np.newaxis,:,np.newaxis])*r[:,i],np.sqrt(m[np.newaxis,:,np.newaxis])*r[:,i]),axis=1)/np.sum(m)
    ev_tmp=[np.linalg.eigh(IT[j])[1][:,-1] for j in range(nfr)] ### eigh function is more suitable for hermitian matrix. The z component of the last column is cos_theta.
    ev[:,i]=ev_tmp*np.sign(np.sum(ev_tmp*drNC[:,i],axis=-1,keepdims=True))
    evr[:,i]=np.sum(rc[:,i]*ev[:,i],axis=-1)/np.linalg.norm(rc[:,i],axis=-1)

bin_r=np.linspace(0,np.linalg.norm(b[0,:,1]-b[0,:,0])/2,nhi_r+1)
bin_mag=np.linspace(-1,1,nhi_mag+1)

mag_mean=np.zeros(nhi_r)
mag_std=np.zeros(nhi_r)
rc=np.linalg.norm(rc,axis=-1)
for i in range(nhi_r):
    idx=(rc>bin_r[i]) & (rc<=bin_r[i+1])
    evr_=evr[idx]
    if len(evr_)>0:
        mag_mean[i]=np.mean(evr_)
        mag_std[i]=np.std(evr_)

mag=np.zeros(nhi_mag)
p=np.zeros((nhi_r,nhi_mag))
hist,xedge,yedge=np.histogram2d(rc.reshape(-1),evr.reshape(-1),bins=[bin_r,bin_mag])
r=(xedge[:-1]+xedge[1:])/2
mag=(yedge[:-1]+yedge[1:])/2
p=hist/np.sum(hist*(np.max(r)-np.min(r))/nhi_r*(np.max(mag)-np.min(mag))/nhi_mag)

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.savez(f'{outname}.npz',r=r,mag=mag,p=p,mag_mean=mag_mean,mag_std=mag_std)