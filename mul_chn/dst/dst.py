"""Import Modules"""
import argparse
import numpy as np
import os
from scipy.spatial.distance import squareform,pdist
from scipy.optimize import curve_fit
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
outname=f'dst/{mdl}_{T}'
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]
nhi=1000

"""Read Data and Calculate"""
m=np.array([m[seq[i]] for i in range(nat)])
r,tp,b=xtc_rd(xtcfile,f'{xtcfile[:-9]}.dat')
z=r[:,:,-1]
nfr,N=np.shape(z)
nch=round(N/nat)
b-=np.mean(b,axis=-1,keepdims=True)
err=1e-5
for i in range(nfr):
    z[i,z[i]<b[i,-1,0]]=(1-err)*b[i,-1,0]+err*b[i,-1,1]
    z[i,z[i]>b[i,-1,1]]=(1-err)*b[i,-1,1]+err*b[i,-1,0]
M=np.tile(m,(nfr,nch))
bin=np.linspace(b[0,-1,0],b[0,-1,1],nhi+1)
Z=(bin[:-1]+bin[1:])/2
bin_idx=np.digitize(z.reshape(-1),bin)-1
hist=np.bincount(bin_idx,weights=M.reshape(-1),minlength=nhi)
dst=hist/nfr/(bin[1]-bin[0])

def rho(z,z1,z2,d,rhol,rhoh):
    return rhol+(rhoh-rhol)/2*(np.tanh((z-z1)/d)-np.tanh((z-z2)/d))
init=[np.min(Z)+0.3*(np.max(Z)-np.min(Z)),np.min(Z)+0.7*(np.max(Z)-np.min(Z)),1.0,np.min(dst),np.max(dst)]
bounds=([np.min(Z),np.min(Z),0.0,0.0,0.0],
        [np.max(Z),np.max(Z),np.inf,np.inf,np.inf])
par,cov=curve_fit(rho,Z,dst,p0=init,bounds=bounds)

"""Output Data"""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
pyw=open(f'{outname}.pyw','w')
pyw.write(f'Z Coordinate, Density:\n')
for i in range(nhi):
    pyw.write(f'{Z[i]} {dst[i]}\n')
pyw.write('\n')
pyw.write(f'Position of Left Interface, Position of Right Interface, Thickness of Interface, Density of Dilute Zone, Density of Dense Zone:\n')
for i in range(5):
    pyw.write(f'{par[i]} ')
pyw.close()