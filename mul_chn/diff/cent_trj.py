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
parser.add_argument('-f1',type=str,help='The trajectory file')
parser.add_argument('-f2',type=str,help='The trajectory file')
args=parser.parse_args()
xtcfile1=args.f1
xtcfile2=args.f2
mdl='_'.join(xtcfile1.split('/')[-2].split('_')[-3:-1])
T=xtcfile1.split('/')[-2].split('_')[-1]
outname=f'cent_trj_81/{mdl}_{T}'
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]
dt=10e-6*10000 ### ns

""""""
r1=xtc_rd(xtcfile1,f'{xtcfile1[:-8]}.dat')[0]
r2=xtc_rd(xtcfile2,f'{xtcfile2[:-12]}.dat')[0]
nfr,N,_=np.shape(r1)
nch=round(N/nat)
r1=r1.reshape(nfr,nch,nat,3)
r2=r2.reshape(nfr,nch,nat,3)
# m=np.array([m[seq[i]] for i in range(nat)])
# M=np.tile(m,(nch,1))
# rc1=np.sum(M[np.newaxis,:,:,np.newaxis]*r1,axis=2)/np.sum(m)
# rc2=np.sum(M[np.newaxis,:,:,np.newaxis]*r2,axis=2)/np.sum(m)
rc1=r1[:,:,81]
rc2=r2[:,:,81]

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.savez(f'{outname}.npz',rc1=rc1,rc2=rc2)