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
outname=f'loc_z/{mdl}_{T}'
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]

"""Read Data"""
m=np.array([m[seq[i]] for i in range(nat)])
r,tp,b=xtc_rd(xtcfile,f'{xtcfile[:-18]}.dat')
z=r[:,:,-1]
nfr,N=np.shape(z)
nch=round(N/nat)

z1,z2,d=pyw(pywfile,'Position')[:3,0]

"""Judge If a Molecule Is inside, in the Surface of, or Outside the Condensate"""
z=z.reshape(nfr,nch,nat)
M=np.tile(m,(nfr,nch,1))
zc=np.sum(M*z,axis=-1)/np.sum(m)

# loc=np.heaviside(zc-z1-d,0.0)-np.heaviside(zc-z2+d,0.0)
loc=np.heaviside(zc-z1-d,0.0)-np.heaviside(zc-z2+d,0.0)+np.heaviside(zc-z1+d,0.0)-np.heaviside(zc-z2-d,0.0)
loc=np.round(loc).astype(int)

# """Create Slice for Each Period of Time Spent inside"""
# slc=[]
# for i in range(nch):
#     diff=loc[1:,i]-loc[:-1,i]
#     if np.sum(diff)==0:
#         diff=np.concatenate((np.array([1]),diff))
#         diff=np.concatenate((diff,np.array([-1])))
#     else:
#         for j in range(nfr-1):
#             if diff[j]!=0:
#                 if diff[j]==-1:
#                     diff=np.concatenate((np.array([1]),diff))
#                 else:
#                     diff=np.concatenate((np.array([0]),diff))
#                 break
#         for j in range(nfr-1,-1,-1):
#             if diff[j]!=0:
#                 if diff[j]==1:
#                     diff=np.concatenate((diff,np.array([-1])))
#                 else:
#                     diff=np.concatenate((diff,np.array([0])))
#                 break
#     slc.append([])
#     for j in range(nfr+1):
#         if diff[j]==1:
#             slc[-1].append([j])
#         if diff[j]==-1:
#             slc[-1][-1]+=[j]

# """Create Increasing Number Identity for Each Period of Time Spent inside"""
# idx_slc=np.zeros((nfr,nch),dtype=int)-1
# for i in range(nch):
#     for j in range(len(slc[i])):
#         idx_slc[slc[i][j][0]:slc[i][j][1],i]=j

"""Create Increasing Number Identity for Each Period of Time Spent Inside, in the Surface, and Outside"""
idx_slc=np.zeros((nfr,nch),dtype=int)-1
for i in range(nch):
    idx_in=0
    idx_surf=0
    idx_out=0
    for j in range(nfr):
        if loc[j,i]==2:
            idx_slc[j,i]=idx_in
            if j<nfr-1:
                if loc[j+1,i]!=2:
                    idx_in+=1
        elif loc[j,i]==1:
            idx_slc[j,i]=idx_surf
            if j<nfr-1:
                if loc[j+1,i]!=1:
                    idx_surf+=1
        elif loc[j,i]==0:
            idx_slc[j,i]=idx_out
            if j<nfr-1:
                if loc[j+1,i]!=0:
                    idx_out+=1

"""Output Data"""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
pyw=open(f'{outname}.pyw','w')
pyw.write(f'Frame Index, Chain Index, Location (in-2 Surface-1 out-0), Increasing Number Identity of Location:\n')
for i in range(nfr):
    for j in range(nch):
        pyw.write(f'{i} {j} {loc[i,j]} {idx_slc[i,j]}\n')
pyw.close()

# import numpy as np
# inout=np.array([0,0,1,2,2,2,0,1,1,2,2,1,1,1,0,0,2])
# nfr=len(inout)
# idx_slc=np.zeros(nfr,dtype=int)-1
# idx_in=0
# idx_surf=0
# idx_out=0
# for j in range(nfr):
#     if inout[j]==2.0:
#         idx_slc[j]=idx_in
#         if j<nfr-1:
#             if inout[j+1]!=2.0:
#                 idx_in+=1
#     elif inout[j]==1.0:
#         idx_slc[j]=idx_surf
#         if j<nfr-1:
#             if inout[j+1]!=1.0:
#                 idx_surf+=1
#     elif inout[j]==0.0:
#         idx_slc[j]=idx_out
#         if j<nfr-1:
#             if inout[j+1]!=0.0:
#                 idx_out+=1

# print(inout)
# print(idx_slc)