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
outname=f'msd/{idx}'
dt=10e-6*10000 ### ns

""""""
r=xtc_rd(xtcfile,f'{xtcfile[:-4]}.dat')[0]
nfr,nat,_=np.shape(r)

d2=np.sum((r[np.newaxis,:,:,:]-r[:,np.newaxis,:,:])**2,axis=-1)
d2mean=np.zeros((nfr,nat))
for i in range(nfr):
    for j in range(nfr-i):
        d2mean[i,:]+=d2[j,i+j,:]
    d2mean[i,:]/=nfr-i

# trim=np.tri(nfr,nfr,k=0) ### Lower triangular matrix.
# d2=d2*trim[:,:,np.newaxis]
# idx=np.abs(np.arange(nfr)-np.arange(nfr)[:,np.newaxis])
# d2sum=np.zeros((nfr,nat))
# np.add.at(d2sum,idx,d2)
# norm=np.arange(nfr,0,-1)
# d2mean=d2sum/norm[:,np.newaxis]

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
pyw=open(f'{outname}.pyw','w')
pyw.write('Lag Time (ns), Residue Index, Mean-squared Displacement:\n')
for i in range(nfr):
    for j in range(nat):
        pyw.write(f'{i*dt} {j} {d2mean[i,j]}\n')