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
parser.add_argument('-f',type=str,help='The .npz file')
args=parser.parse_args()
npzfile=args.f
mdl='_'.join(npzfile.split('/')[-1].split('_')[-3:-1])
T=float(npzfile.split('/')[-1].split('_')[-1][:-4])
outname=f'grf_mod/{mdl}_{T}'
lx=150.0
ly=150.0
lz=500.0
kB=1.987e-3
dt=10 ### fs
unit=1.458e-5
# Pcut=50.0

""""""
P=np.load(npzfile)['Pchn']
# P=np.heaviside(P-Pcut,0.0)
deg=np.sum(P,axis=1)
nfr,nch=np.shape(P)[:2]
tau=np.zeros((nfr,nch-1))
for i in range(nfr):
    Z=-P[i]
    np.fill_diagonal(Z,deg[i])
    eval=np.linalg.eigvalsh(Z)[1:] ### eigvalsh for symmetric matrix
    tau[i]=1/eval/T

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.save(f'{outname}.npy',tau)
                    