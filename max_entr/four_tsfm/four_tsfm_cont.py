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
outname=f'four_tsfm_cont/{idx}'
mu=1.0
csgm=float(idx.split('_')[0][:-3])
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
sgm=[5.04,6.56,5.68,5.58,5.48,6.02,5.92,4.50,6.08,6.18,6.18,6.36,6.18,6.36,5.56,5.18,5.62,6.78,6.46,5.86]
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]
dt=10e-6*10000 ### ns

""""""
r=xtc_rd(xtcfile,f'{xtcfile[:-4]}.dat')[0]
nfr,nat,_=np.shape(r)
D=np.zeros((nfr,nat,nat))
for i in range(nfr):
    D[i]=squareform(pdist(r[i],'euclidean'))

Sgm=np.zeros(nat)
for i in range(nat):
    Sgm[i]=sgm[seq[i]]
Sgm=(Sgm[np.newaxis,:]+Sgm[:,np.newaxis])/2
Q=0.5*(1-np.tanh(mu*(D-(csgm*Sgm)[np.newaxis])))

F=np.zeros((nfr,nat,nat),dtype=complex)
for i in range(nat):
    for j in range(i):
        F[:,i,j]=np.fft.fft(Q[:,i,j])
        F[:,j,i]=F[:,i,j]
freq=np.fft.fftfreq(nfr,dt)
A=np.abs(F)
phi=np.angle(F)
freq=freq[:round(nfr/2)]
A=A[:round(nfr/2)]
phi=phi[:round(nfr/2)]

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.savez(f'{outname}.npz',freq=freq,A=A,phi=phi)