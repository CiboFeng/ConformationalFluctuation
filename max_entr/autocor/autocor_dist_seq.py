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
outname=f'autocor_dist_seq/{idx}'
dt=10e-6*10000 ### ns

""""""
r=xtc_rd(xtcfile,f'{xtcfile[:-4]}.dat')[0]
nfr,nat,_=np.shape(r)
D=np.zeros((nfr,nat,nat))
for i in range(nfr):
    D[i]=squareform(pdist(r[i],'euclidean'))

d=np.zeros((nfr,nat))
for i in range(nat):
    for j in range(nat-i):
        d[:,i]+=D[:,j,i+j]
    d[:,i]/=nat-i

dmean=np.mean(d,axis=0,keepdims=True)
dvar=np.var(d,axis=0,keepdims=True)
dc=d-dmean
dcrs=dc[np.newaxis,:,:]*dc[:,np.newaxis,:]
dcrsmean=np.zeros((nfr,nat))
for i in range(nfr):
    for j in range(nfr-i):
        dcrsmean[i]+=dcrs[j,i+j]
    dcrsmean[i]/=nfr-i
ac=dcrsmean/dvar

# trim=np.tri(nfr,nfr,k=0) ### Lower triangular matrix.
# dcrs=dcrs*trim[:,:,np.newaxis]
# idx=np.abs(np.arange(nfr)-np.arange(nfr)[:,np.newaxis])
# dcrssum=np.zeros(nfr,nat)
# np.add.at(dcrssum,idx,dcrs)
# norm=np.arange(nfr,0,-1)
# dcrsmean=dcrssum/norm[:,np.newaxis]
# ac=dcrsmean/dvar

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
pyw=open(f'{outname}.pyw','w')
pyw.write('Lag Time (ns), Sequence Distance, Autocorrelation of End-end Distance:\n')
for i in range(nfr):
    for j in range(nat):
        pyw.write(f'{i*dt} {j} {ac[i,j]}\n')