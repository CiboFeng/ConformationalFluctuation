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
idx='_'.join(xtcfile.split('/')[-2].split('_')[-3:])
outname=f'autocor_dist/{idx}'
dt=10e-6*10000 ### ns

""""""
r=xtc_rd(xtcfile,f'{xtcfile[:-4]}.dat')[0]
nfr,nat,_=np.shape(r)
D=np.zeros((nfr,nat,nat))
for i in range(nfr):
    D[i]=squareform(pdist(r[i],'euclidean'))

# Dmean=np.mean(D,axis=0,keepdims=True)
# Dvar=np.var(D,axis=0,keepdims=True)
# Dc=D-Dmean
# Dcrs=Dc[np.newaxis,:,:,:]*Dc[:,np.newaxis,:,:]
# Dcrsmean=np.zeros((nfr,nat,nat))
# for i in range(nfr):
#     for j in range(nfr-i):
#         Dcrsmean[i,:,:]+=Dcrs[j,i+j,:,:]
#     Dcrsmean[i,:,:]/=nfr-i
# ac=Dcrsmean/Dvar

# ac=np.ones((nfr,nat,nat))
# for i in range(nat):
#     for j in range(i):
#         Dmean=np.mean(D[:,i,j])
#         Dvar=np.var(D[:,i,j])
#         Dc=D[:,i,j]-Dmean
#         Dcrs=Dc[np.newaxis,:]*Dc[:,np.newaxis]
#         Dcrsmean=np.zeros(nfr)
#         for k in range(nfr):
#             for l in range(nfr-k):
#                 Dcrsmean[k]+=Dcrs[l,k+l]
#             Dcrsmean[k]/=nfr-k
#         ac[:,i,j]=Dcrsmean/Dvar
#         ac[:,j,i]=ac[:,i,j]

Dmean=np.mean(D,axis=0,keepdims=True)
Dvar=np.var(D,axis=0,keepdims=True)
Dc=D-Dmean
Dcrsmean=np.zeros((nfr,nat,nat))
for i in range(nfr):
    for j in range(nfr-i):
        Dcrsmean[i]+=Dc[j]*Dc[i+j]
    Dcrsmean[i]/=nfr-i
ac=Dcrsmean/Dvar

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.save(f'{outname}.npy',ac)