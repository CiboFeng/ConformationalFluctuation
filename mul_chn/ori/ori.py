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
outname=f'ori/{mdl}_{T}'
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]
nhi_z=250
nhi_tht=50

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
tht=np.zeros((nfr,nch))
thtabs=np.zeros((nfr,nch))
for i in range(nch):
    IT=np.sum(np.einsum('ijk,ijl->ijkl',np.sqrt(m[np.newaxis,:,np.newaxis])*r[:,i],np.sqrt(m[np.newaxis,:,np.newaxis])*r[:,i]),axis=1)/np.sum(m)
    ev_tmp=np.array([np.linalg.eigh(IT[j])[1][:,-1] for j in range(nfr)]) ### eigh function is more suitable for hermitian matrix. The z component of the last column is cos_theta.
    evz=ev_tmp[:,-1]
    evz=evz*np.sign(np.sum(ev_tmp*drNC[:,i],axis=-1))
    tht[:,i]=np.arccos(np.clip(evz,-1.0,1.0))
    thtabs[:,i]=np.arccos(np.clip(np.abs(evz),0.0,1.0))

bin_z=np.linspace((b[0,-1,0]-b[0,-1,1])/2,(b[0,-1,1]-b[0,-1,0])/2,nhi_z+1)
bin_tht=np.linspace(0,np.pi,nhi_tht+1)
bin_thtabs=np.linspace(0,np.pi/2,nhi_tht+1)

tht_mean=np.zeros(nhi_z)
tht_std=np.zeros(nhi_z)
thtabs_mean=np.zeros(nhi_z)
thtabs_std=np.zeros(nhi_z)
for i in range(nhi_z):
    idx=(zc>bin_z[i]) & (zc<=bin_z[i+1]) ### Can't be tht[binz[i]<zc<=binz[i+1]]
    tht_=tht[idx]
    thtabs_=thtabs[idx]
    if len(tht_)>0:
        tht_mean[i]=np.mean(tht_)
        tht_std[i]=np.std(tht_)
        thtabs_mean[i]=np.mean(thtabs_)
        thtabs_std[i]=np.std(thtabs_)

hist,xedge,yedge=np.histogram2d(zc.reshape(-1),tht.reshape(-1),bins=[bin_z,bin_tht])
z=(xedge[:-1]+xedge[1:])/2
tht=(yedge[:-1]+yedge[1:])/2
p=hist/np.sum(hist*(np.max(z)-np.min(z))/nhi_z*(np.max(tht)-np.min(tht))/nhi_tht)
hist,xedge,yedge=np.histogram2d(zc.reshape(-1),thtabs.reshape(-1),bins=[bin_z,bin_thtabs])
thtabs=(yedge[:-1]+yedge[1:])/2
pabs=hist/np.sum(hist*(np.max(z)-np.min(z))/nhi_z*(np.max(thtabs)-np.min(thtabs))/nhi_tht)

# z1,z2,d=pyw(pywfile,'Position')[:3,0]
# binz=np.linspace(z1-d,z1+d,nhiz+1)

# Tht=np.zeros((nhiz,nhitht))
# p_Tht=np.zeros((nhiz,nhitht))
# Tht_mean=np.zeros(nhiz)
# Tht_std=np.zeros(nhiz)
# bin=np.linspace(0,np.pi/2,nhitht+1)
# for i in range(nhiz):
#     tht_=tht[binz[i]<zc<=binz[i+1]]
#     hist,edge=np.histogram(tht_,bins=bin)
#     Tht[i]=(edge[:-1]+edge[1:])/2
#     p_Tht[i]=hist/np.sum(hist*np.pi/2/nhitht)
#     Tht_mean[i]=np.mean(tht_)
#     Tht_std[i]=np.std(tht_)

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.savez(f'{outname}.npz',z=z,tht=tht,p=p,tht_mean=tht_mean,tht_std=tht_std,thtabs=thtabs,pabs=pabs,thtabs_mean=thtabs_mean,thtabs_std=thtabs_std)