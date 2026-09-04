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
outname=f'cont_seq_life/{idx}'
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
nhi=20

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
# Q=0.5*(1-np.tanh(mu*(D-(csgm*Sgm)[np.newaxis])))
Q=1-np.heaviside(D-(csgm*Sgm)[np.newaxis],0.0)

Qpad=np.pad(Q,((1,1),(0,0),(0,0)),'constant',constant_values=-1)
dQ=Qpad[1:]-Qpad[:-1]
Qs=np.where(dQ[:-1]!=0,1,0)
Qe=np.where(dQ[1:]!=0,1,0)

bin=np.linspace(-1,3,nhi+1)
T0=np.zeros((nat,nhi))
p_T0=np.zeros((nat,nhi))
T0_mean=np.zeros(nat)
T0_std=np.zeros(nat)
T1=np.zeros((nat,nhi))
p_T1=np.zeros((nat,nhi))
T1_mean=np.zeros(nat)
T1_std=np.zeros(nat)
for i in range(nat):
    seglen=[]
    segval=[]
    for j in range(nat-i):
        idxs=np.where(Qs[:,j,i+j]==1)[0]
        idxe=np.where(Qe[:,j,i+j]==1)[0]
        for k,l in zip(idxs,idxe):
            seglen+=[l-k+1]
            segval+=[Q[k,j,i+j]]
    seglen=np.array(seglen)
    segval=np.array(segval)
    len0=seglen[segval==0]
    len1=seglen[segval==1]
    len0=np.log10(np.array([x for x in len0 if x!=0])*dt)
    len1=np.log10(np.array([x for x in len1 if x!=0])*dt)
    hist,edge=np.histogram(len0,bins=bin)
    # hist,edge=np.histogram(len0,bins=nhi)
    T0[i]=(edge[:-1]+edge[1:])/2
    p_T0[i]=hist/np.sum(hist*(np.max(T0[i])-np.min(T0[i]))/nhi)
    T0_mean[i]=np.mean(len0)
    T0_std[i]=np.std(len0)
    hist,edge=np.histogram(len1,bins=bin)
    # hist,edge=np.histogram(len1,bins=nhi)
    T1[i]=(edge[:-1]+edge[1:])/2
    p_T1[i]=hist/np.sum(hist*(np.max(T1[i])-np.min(T1[i]))/nhi)
    T1_mean[i]=np.mean(len1)
    T1_std[i]=np.std(len1)

# q=np.zeros((nfr,nat))
# for i in range(nat):
#     for j in range(nat-i):
#         q[:,i]+=Q[:,j,i+j]
#     q[:,i]/=nat-i
# q=np.heaviside(q-0.5,0.0)

# qpad=np.pad(q,((1,1),(0,0)),'constant',constant_values=-1)
# dq=qpad[1:]-qpad[:-1]
# qs=np.where(dq[:-1]!=0,1,0)
# qe=np.where(dq[1:]!=0,1,0)

# bin=np.linspace(-1,3,nhi+1)
# T0=np.zeros((nat,nhi))
# p_T0=np.zeros((nat,nhi))
# T0_mean=np.zeros(nat)
# T0_std=np.zeros(nat)
# T1=np.zeros((nat,nhi))
# p_T1=np.zeros((nat,nhi))
# T1_mean=np.zeros(nat)
# T1_std=np.zeros(nat)
# for i in range(nat):
#     seglen=[]
#     segval=[]
#     idxs=np.where(qs[:,i]==1)[0]
#     idxe=np.where(qe[:,i]==1)[0]
#     for j,k in zip(idxs,idxe):
#         seglen+=[k-j+1]
#         segval+=[q[j,i]]
#     seglen=np.array(seglen)
#     segval=np.array(segval)
#     len0=seglen[segval==0]
#     len1=seglen[segval==1]
#     len0=np.log10(np.array([x for x in len0 if x!=0])*dt)
#     len1=np.log10(np.array([x for x in len1 if x!=0])*dt)
#     hist,edge=np.histogram(len0,bins=bin)
#     # hist,edge=np.histogram(len0,bins=nhi)
#     T0[i]=(edge[:-1]+edge[1:])/2
#     p_T0[i]=hist/np.sum(hist*(np.max(T0[i])-np.min(T0[i]))/nhi)
#     T0_mean[i]=np.mean(len0)
#     T0_std[i]=np.std(len0)
#     hist,edge=np.histogram(len1,bins=bin)
#     # hist,edge=np.histogram(len1,bins=nhi)
#     T1[i]=(edge[:-1]+edge[1:])/2
#     p_T1[i]=hist/np.sum(hist*(np.max(T1[i])-np.min(T1[i]))/nhi)
#     T1_mean[i]=np.mean(len1)
#     T1_std[i]=np.std(len1)
        
""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.savez(f'{outname}.npz',T0=T0,p_T0=p_T0,T0_mean=T0_mean,T0_std=T0_std,
         T1=T1,p_T1=p_T1,T1_mean=T1_mean,T1_std=T1_std)