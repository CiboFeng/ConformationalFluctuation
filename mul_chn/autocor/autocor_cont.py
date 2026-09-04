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
pywdir='../loc_z/loc_z'
mdl='_'.join(xtcfile.split('/')[-2].split('_')[-3:-1])
T=xtcfile.split('/')[-2].split('_')[-1]
outname=f'autocor_dist_in/{mdl}_{T}'
nat=163
dt=10e-6*10000 ### ns

""""""
r,tp,b=xtc_rd(xtcfile,f'{xtcfile[:-14]}.dat')
nfr,N,_=np.shape(r)
nch=round(N/nat)
r=r.reshape(nfr,nch,nat,3)
D=np.zeros((nfr,nch,nat,nat))
for i in range(nfr):
    for j in range(nch):
        D[i,j]=squareform(pdist(r[i,j],'euclidean'))
Dmean=np.mean(D,axis=(0,1),keepdims=True)
Dvar=np.var(D,axis=(0,1))[np.newaxis]
Dc=D-Dmean
Dcrsmean=np.zeros((nfr,nat,nat))

loc,idx_slc=pyw(f'{pywdir}/{mdl}_{T}.pyw','Location',int)[2:]
loc=loc.reshape(nfr,nch)
idx_slc=idx_slc.reshape(nfr,nch)

for i in range(nfr):
    num=0
    for j in range(nfr-i):
        for k in range(nch):
            if loc[j,k]==2 and loc[i+j,k]==2 and idx_slc[j,k]==idx_slc[i+j,k]:
                Dcrsmean[i]+=Dc[j,k]*Dc[i+j,k]
                num+=1
    if num>0:
        Dcrsmean[i]/=num
ac=Dcrsmean/Dvar

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.save(f'{outname}.npy',ac)