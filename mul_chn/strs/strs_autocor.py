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
parser.add_argument('-f',type=str,help='The log file')
args=parser.parse_args()
logfile=args.f
mdl='_'.join(logfile.split('/')[-2].split('_')[-3:-1])
T=float(logfile.split('/')[-2].split('_')[-1])
outname=f'strs_autocor/{mdl}_{T}'
lx=150.0
ly=150.0
lz=500.0
kB=1.987e-3
dt=10 ### fs
unit=1.458e-5

""""""
P=[[] for _ in range(6)]
log=open(logfile)
lslog=log.readlines()
for i in range(len(lslog)):
    if lslog[i][:16]=='Step Pxx Pyy Pzz':
        for j in range(i+1,len(lslog)):
            ls_log=lslog[j].strip('\n').split()
            try:
                float(ls_log[0])
            except Exception:
                break
            for k in range(6):
                P[k]+=[float(ls_log[k+1])]
P=np.array(P)*unit
nfr=np.shape(P)[1]

Pcrs=np.zeros((6,nfr))
for i in range(nfr):
    for j in range(nfr-i):
        Pcrs[:,i]+=P[:,j]*P[:,i+j]
    Pcrs[:,i]/=nfr-i
ac=Pcrs*lx*ly*lz/(kB*T)

nfr=1000
ac=ac[:,:nfr]

G=np.zeros((6,nfr),dtype=complex)
for i in range(6):
    G[i]=np.fft.fft(ac[i])*dt
freq=np.fft.fftfreq(nfr,dt)
A=np.abs(G)
phi=np.angle(G)
freq=freq[:round(nfr/2)]
A=A[:,:round(nfr/2)]
phi=phi[:,:round(nfr/2)]

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.savez(f'{outname}.npz',ac=ac,freq=freq,A=A,phi=phi)
                    