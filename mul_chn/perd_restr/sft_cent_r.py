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
resol=0.1
S=0.7

"""Read Data"""
m=np.array([m[seq[i]] for i in range(nat)])
r,tp,b=xtc_rd(ixtcfile,f'{ixtcfile[:-8]}.dat')
nfr,N,_=np.shape(r)
nfr=round(nfr*S)
r=r[-nfr:]
b=b[-nfr:]
nch=round(N/nat)
M=np.tile(m,nch)
b-=np.mean(b,axis=-1,keepdims=True)

"""Shift the Condensate to Center of Box According to Center of Mass"""
rc=np.sum(M[np.newaxis,:,np.newaxis]*r,axis=1)/np.sum(M)
r=(sft(r,b+rc[:,:,np.newaxis],[[0,N]])-rc[:,np.newaxis,:])

rcs=[rc]
idx=list(range(nfr))
i=1
while len(idx)>0:
    rc=np.sum(M[np.newaxis,:,np.newaxis]*r,axis=1)/np.sum(M)
    rcs+=[rc]
    idx=np.where(np.linalg.norm(rcs[i]-rcs[i-1],axis=-1)>resol)[0].tolist()
    for j in idx:
        print(i,j)
        r[j:j+1]=(sft(r[j:j+1],b[j:j+1]+rc[j:j+1,:,np.newaxis],[[0,N]])-rc[j:j+1,np.newaxis,:])
    i+=1

"""Output Data"""
xtc_wrt(oxtcfile,r,tp,b)
