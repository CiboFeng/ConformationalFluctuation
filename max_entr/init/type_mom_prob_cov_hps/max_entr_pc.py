"""Import Modules"""
import argparse
import numpy as np
from scipy.spatial.distance import squareform,pdist
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/opt/mypylib')
sys.path.append('/hpc2hdd/home/chuo657/cbfengphy/opt/mypylib')
from pyw import pyw
from trj import trj_rd

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-s',type=str,help='The sequence file')
parser.add_argument('-t',type=str,help='The .trj file')
parser.add_argument('-p',type=str,help='The .npz file containing contact probabilities and covariance matrix')
parser.add_argument('-n',type=int,help='The number of moment')
args=parser.parse_args()
seqfile=args.s
trjfile=args.t
pcfile=args.p
nmo=args.n
sgm=[5.04,6.56,5.68,5.58,5.48,6.02,5.92,4.50,6.08,6.18,6.18,6.36,6.18,6.36,5.56,5.18,5.62,6.78,6.46,5.86]
mu=1.0
csgm=1.5
ntp=20

"""Read Sequence File and .trj File, and Calculate"""
seq=pyw(seqfile,'Sequence','int')[0]
N=len(seq)
r=trj_rd(trjfile)[0]
nfr=np.shape(r)[0]

D=np.zeros((nfr,N,N))
for i in range(nfr):
    D[i]=squareform(pdist(r[i],'euclidean'))
Sgm=np.zeros(N)
for i in range(N):
    Sgm[i]=sgm[seq[i]]
Sgm=(Sgm[np.newaxis,:]+Sgm[:,np.newaxis])/2
P0=(0.5*(1-np.tanh(mu*(D-csgm*Sgm[np.newaxis,:,:]))))
P=np.zeros((nfr,N,N,nmo))
for i in range(nmo):
    P[:,:,:,i]=P0**(i+1)
for i in range(N):
    for j in range(N):
        if abs(i-j)<4:
            P[:,i,j,:]=0
p=np.zeros((nfr,ntp,ntp,nmo))
for i in range(N):
    for j in range(i+1):
        p[:,seq[i],seq[j],:]+=P[:,i,j,:]
        if seq[i]!=seq[j]:
            p[:,seq[j],seq[i],:]+=P[:,i,j,:]
idx=np.tril_indices(ntp)
p=p[:,idx[0],idx[1]]
p=p.transpose(0,2,1).reshape(nfr,-1)
C=np.dot(p.T,p)/nfr
p=np.mean(p,axis=0)
P=np.mean(P,axis=0)

"""Write to .npz File"""
np.savez(pcfile,nfr=nfr,p=p,C=C,P=P)