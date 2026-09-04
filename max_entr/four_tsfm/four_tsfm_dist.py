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
outname=f'four_tsfm_dist/{idx}'
dt=10e-6*10000 ### ns

""""""
r=xtc_rd(xtcfile,f'{xtcfile[:-4]}.dat')[0]
nfr,nat,_=np.shape(r)
D=np.zeros((nfr,nat,nat))
for i in range(nfr):
    D[i]=squareform(pdist(r[i],'euclidean'))

F=np.zeros((nfr,nat,nat),dtype=complex)
for i in range(nat):
    for j in range(i):
        F[:,i,j]=np.fft.fft(D[:,i,j])
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