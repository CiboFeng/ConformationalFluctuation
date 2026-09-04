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
parser.add_argument('-f',type=str,help='The trajectory file')
args=parser.parse_args()
xtcfile=args.f
# npzdir='../loc_z/loc_z_res'
pywdir='../loc_z/loc_z'
mdl='_'.join(xtcfile.split('/')[-2].split('_')[-3:-1])
T=xtcfile.split('/')[-2].split('_')[-1]
outname=f'msd_res_in/{mdl}_{T}'
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]

""""""
r=xtc_rd(xtcfile,f'{xtcfile[:-12]}.dat')[0]
nfr,N,_=np.shape(r)
nch=round(N/nat)
r=r.reshape(nfr,nch,nat,3)
m=np.array([m[seq[i]] for i in range(nat)])
Rc=np.sum(m[np.newaxis,np.newaxis,:,np.newaxis]*r,axis=(1,2))/np.sum(m)/nch
r[:,:,:,-1]-=Rc[:,np.newaxis,np.newaxis,-1]

# q=np.load(f'{npzdir}/{idx}_{T}.npz')
# loc=q['loc']
# idx_slc=q['idx_slc']
loc,idx_slc=pyw(f'{pywdir}/{mdl}_{T}.pyw','Location',int)[2:]
loc=loc.reshape(nfr,nch)
idx_slc=idx_slc.reshape(nfr,nch)

d2mean=np.zeros((nfr,nch,nat,3))
# for i in range(nch):
#     for j in range(nat):
#         for k in range(nfr):
#             num=0
#             for l in range(nfr-k):
#                 if loc[l,i,j]==2 and loc[k+l,i,j]==2 and idx_slc[l,i,j]==idx_slc[k+l,i,j]:
#                     d2mean[k,i,j]+=np.sum((r[l,i,j]-r[k+l,i,j])**2)
#                     num+=1
#             if num>0:
#                 d2mean[k,i,j]/=num
for i in range(nch):
    for j in range(nfr):
        num=0
        for k in range(nfr-j):
            if loc[k,i]==2 and loc[j+k,i]==2 and idx_slc[k,i]==idx_slc[j+k,i]:
                d2mean[j,i]+=(r[k,i]-r[j+k,i])**2
                num+=1
        if num>0:
            d2mean[j,i]/=num

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.save(f'{outname}.npy',d2mean)