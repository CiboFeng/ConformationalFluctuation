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
bin_mag=np.linspace(-1,1,nhi_mag+1)

mag_mean=np.zeros((nhi_z,3))
mag_std=np.zeros((nhi_z,3))
for i in range(nhi_z):
    idx=(zc>bin_z[i]) & (zc<=bin_z[i+1]) ### Can't be tht[binz[i]<zc<=binz[i+1]]
    ev_=ev[idx]
    if len(ev_)>0:
        mag_mean[i]=np.mean(ev_,axis=0)
        mag_std[i]=np.std(ev_,axis=0)

mag=np.zeros((nhi_mag,3))
p=np.zeros((nhi_z,nhi_mag,3))
for i in range(3):
    hist,xedge,yedge=np.histogram2d(zc.reshape(-1),ev[:,:,i].reshape(-1),bins=[bin_z,bin_mag])
    z=(xedge[:-1]+xedge[1:])/2
    mag[:,i]=(yedge[:-1]+yedge[1:])/2
    p[:,:,i]=hist/np.sum(hist*(np.max(z)-np.min(z))/nhi_z*(np.max(mag[:,i])-np.min(mag[:,i]))/nhi_mag)

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.savez(f'{outname}.npz',z=z,mag=mag,p=p,mag_mean=mag_mean,mag_std=mag_std)