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
from xtc import xtc_rd,xtc_wrt
from perd_restr import sft

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-i',type=str,help='The inputting trajectory file')
parser.add_argument('-o',type=str,help='The outputting trajectory file')
args=parser.parse_args()
ixtcfile=args.i
oxtcfile=args.o
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]
nhi=1000

"""Read Data"""
m=np.array([m[seq[i]] for i in range(nat)])
r,tp,b=xtc_rd(ixtcfile,f'{ixtcfile[:-8]}.dat')
z=r[:,:,-1]
nfr,N=np.shape(z)
nch=round(N/nat)
M=np.tile(m,nch)
b-=np.mean(b,axis=-1,keepdims=True)

"""Preliminerily Shift the Condensate to Center of Box"""
bin=np.linspace(b[0,-1,0],b[0,-1,1],round(nhi/10)+1) ### linspace contains the end.
for i in range(nfr):
    z[i,z[i]<b[i,-1,0]]=b[i,-1,0]
    z[i,z[i]>b[i,-1,1]]=b[i,-1,1]
    bin_idx=np.digitize(z[i],bin)-1
    hist=np.bincount(bin_idx,weights=M,minlength=round(nhi/10))
    z0_idx=np.where(hist==np.max(hist))[0]
    z0=(bin[z0_idx]+bin[z0_idx+1])/2
    z[i:i+1]=(sft(z[i:i+1,:,np.newaxis],b[i:i+1,2:3]+z0,[[0,N]])-z0)[:,:,0]

"""Precisly Shift the Condensate to Center of Box According to Center of Mass"""
zc=np.sum(M[np.newaxis]*z,axis=1)/np.sum(M)
z=(sft(z[:,:,np.newaxis],b[:,2:3]+zc[:,np.newaxis,np.newaxis],[[0,N]])-zc[:,np.newaxis,np.newaxis])[:,:,0] ### Here a new array z is created. It is not related to r again.
r[:,:,-1]=z

"""Output Data"""
xtc_wrt(oxtcfile,r,tp,b)
