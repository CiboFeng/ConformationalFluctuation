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
mdl='_'.join(xtcfile.split('/')[-2].split('_')[-2:])
outname=f'dst_r/{mdl}'
nrep=20
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
r=[]
b=[]
for i in range(nrep):
    xtcfile_=xtcfile.replace('__',f'_{i}_')
    r_,tp,b_=xtc_rd(xtcfile_,f'{xtcfile_[:-9]}.dat')
    r+=[r_]
    b+=[b_]
r=np.linalg.norm(np.concatenate(r,axis=0),axis=-1)
b=np.concatenate(b,axis=0)
nfr,N=np.shape(r)
nch=round(N/nat)
b-=np.mean(b,axis=-1,keepdims=True)
M=np.tile(m,(nfr,nch))
bin=np.linspace(0,np.linalg.norm(b[0,:,1]-b[0,:,0])/2,nhi+1)
R=(bin[:-1]+bin[1:])/2
bin_idx=np.digitize(r.reshape(-1),bin)-1
hist=np.bincount(bin_idx,weights=M.reshape(-1),minlength=nhi)
dst=hist/nfr/(bin[1]-bin[0])

def rho(r,rs,d,rhol,rhoh):
    return (rhol+(rhoh-rhol)/2*(1-np.tanh((r-rs)/d)))*4*np.pi*r**2
init=[np.min(R)+0.5*(np.max(R)-np.min(R)),1.0,np.min(dst),np.max(dst)]
bounds=([np.min(R),0.0,0.0,0.0],
        [np.max(R),np.inf,np.inf,np.inf])
par,cov=curve_fit(rho,R,dst,p0=init,bounds=bounds)

"""Output Data"""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
pyw=open(f'{outname}.pyw','w')
pyw.write(f'Radial Position, Density:\n')
for i in range(nhi):
    pyw.write(f'{R[i]} {dst[i]/(4*np.pi*R[i]**2)}\n')
pyw.write('\n')
pyw.write(f'Position of Interface, Thickness of Interface, Density of Dilute Zone, Density of Dense Zone:\n')
for i in range(4):
    pyw.write(f'{par[i]} ')
pyw.close()