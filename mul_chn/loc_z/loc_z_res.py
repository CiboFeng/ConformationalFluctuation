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
from pyw import pyw
from xtc import xtc_rd
from perd_restr import sft,str_restr

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-f',type=str,help='The trajectory file')
args=parser.parse_args()
xtcfile=args.f
mdl='_'.join(xtcfile.split('/')[-2].split('_')[-3:-1])
T=xtcfile.split('/')[-2].split('_')[-1]
pywfile=f'../dst/dst/{mdl}_{T}.pyw'
outname=f'loc_z_res/{mdl}_{T}'
nat=163

"""Read Data"""
r,tp,b=xtc_rd(xtcfile,f'{xtcfile[:-18]}.dat')
z=r[:,:,-1]
nfr,N=np.shape(z)
nch=round(N/nat)

z1,z2,d=pyw(pywfile,'Position')[:3,0]

"""Judge If a Residue Is inside, in the Surface of, or Outside the Condensate"""
z=z.reshape(nfr,nch,nat)
loc=np.heaviside(z-z1-d,0.0)-np.heaviside(z-z2+d,0.0)+np.heaviside(z-z1+d,0.0)-np.heaviside(z-z2-d,0.0)
loc=np.round(loc).astype(int)

"""Create Increasing Number Identity for Each Period of Time Spent Inside, in the Surface, and Outside"""
idx_slc=np.zeros((nfr,nch,nat),dtype=int)-1
for i in range(nch):
    for j in range(nat):
        idx_in=0
        idx_surf=0
        idx_out=0
        for k in range(nfr):
            if loc[k,i,j]==2:
                idx_slc[k,i,j]=idx_in
                if k<nfr-1:
                    if loc[k+1,i,j]!=2:
                        idx_in+=1
            elif loc[k,i,j]==1:
                idx_slc[k,i,j]=idx_surf
                if k<nfr-1:
                    if loc[k+1,i,j]!=1:
                        idx_surf+=1
            elif loc[k,i,j]==0:
                idx_slc[k,i,j]=idx_out
                if k<nfr-1:
                    if loc[k+1,i,j]!=0:
                        idx_out+=1

"""Output Data"""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.savez(f'{outname}.npz',loc=loc,idx_slc=idx_slc)