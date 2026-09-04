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
parser.add_argument('-f1',type=str,help='The trajectory file')
parser.add_argument('-f2',type=str,help='The trajectory file')
args=parser.parse_args()
xtcfile1=args.f1
xtcfile2=args.f2
mdl='_'.join(xtcfile1.split('/')[-2].split('_')[-3:-1])
T=xtcfile1.split('/')[-2].split('_')[-1]
outname=f'z_distr/{mdl}_{T}'
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
r,tp,b=xtc_rd(xtcfile1,f'{xtcfile1[:-9]}.dat')
b-=np.mean(b,axis=-1,keepdims=True)

bin=np.linspace(b[0,-1,0],b[0,-1,1],nhi+1)
hist,edge=np.histogram(r[:,:,-1],bins=bin)
z=(edge[:-1]+edge[1:])/2
pz=hist/np.sum(hist*(np.max(z)-np.min(z))/nhi)

r,tp,b=xtc_rd(xtcfile2,f'{xtcfile2[:-18]}.dat')
nfr,N,_=np.shape(r)
nch=round(N/nat)
b-=np.mean(b,axis=-1,keepdims=True)
r=r.reshape(nfr,nch,nat,3)
zc=np.sum(m[np.newaxis,np.newaxis]*r[:,:,:,-1],axis=-1)/np.sum(m)

hist,edge=np.histogram(zc,bins=bin)
zc=(edge[:-1]+edge[1:])/2
pzc=hist/np.sum(hist*(np.max(zc)-np.min(zc))/nhi)

"""Output Data"""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
pyw=open(f'{outname}.pyw','w')
pyw.write(f'Z Coordinate of All Residues, Its Probability Density, Z Coordinate of Mass Center of Each Protein, Its Probability Density:\n')
for i in range(nhi):
    pyw.write(f'{z[i]} {pz[i]} {zc[i]} {pzc[i]}\n')
pyw.close()