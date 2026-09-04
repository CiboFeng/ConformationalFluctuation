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
parser.add_argument('-f',type=str,help='The trajectory file')
args=parser.parse_args()
xtcfile=args.f
mdl='_'.join(xtcfile.split('/')[-2].split('_')[-3:-1])
T=xtcfile.split('/')[-2].split('_')[-1]
outname=f'cont_slab/{mdl}_{T}'
mu=1.0
csgm=float(mdl.split('_')[0][:-3])
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
sgm=[5.04,6.56,5.68,5.58,5.48,6.02,5.92,4.50,6.08,6.18,6.18,6.36,6.18,6.36,5.56,5.18,5.62,6.78,6.46,5.86]
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]

""""""
r,tp,b=xtc_rd(xtcfile,f'{xtcfile[:-8]}.dat')
nfr,N,_=np.shape(r)
nch=round(N/nat)
Sgm=np.zeros(nat)
for i in range(nat):
    Sgm[i]=sgm[seq[i]]
Sgm=np.tile(Sgm,nch)
Sgm=(Sgm[np.newaxis,:]+Sgm[:,np.newaxis])/2

Pintra=np.zeros((nat,nat))
Pinter=np.zeros((nat,nat))
Pchn=np.zeros((nfr,nch,nch))
for i in range(nfr):
    # p=0.5*(1-np.tanh(mu*(squareform(pdist(r[i],'euclidean'))-csgm*Sgm)))
    dx=r[i,:,np.newaxis,0]-r[i,np.newaxis,:,0] ### Dividing into xyz can accelerate.
    dy=r[i,:,np.newaxis,1]-r[i,np.newaxis,:,1]
    dz=r[i,:,np.newaxis,2]-r[i,np.newaxis,:,2]
    lx=b[i,0,1]-b[i,0,0]
    ly=b[i,1,1]-b[i,1,0]
    lz=b[i,2,1]-b[i,2,0]
    dx-=np.round(dx/lx)*lx
    dy-=np.round(dy/ly)*ly
    dz-=np.round(dz/lz)*lz
    p=0.5*(1-np.tanh(mu*(np.sqrt(dx**2+dy**2+dz**2)-csgm*Sgm))) ### np.sqrt(dx**2+dy**2+dz**2) is faster than np.sqrt(np.sum(d**2,axis=-1)).
    p=p.reshape(nch,nat,nch,nat)
    pintra=np.zeros((nat,nat))
    pinter=np.zeros((nat,nat))
    for j in range(nch):
        pintra+=p[j,:,j,:]
        for k in range(j):
            pinter+=p[j,:,k,:]
    Pintra+=pintra/nch
    Pinter+=pinter/(nch*(nch-1)/2)
    Pchn[i]=np.sum(p,axis=(1,3))
Pintra/=nfr
Pinter/=nfr

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.savez(f'{outname}.npz',Pintra=Pintra,Pinter=Pinter,Pchn=Pchn)
