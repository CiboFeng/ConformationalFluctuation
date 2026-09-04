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
parser.add_argument('-f',type=str,help='The .npz file')
args=parser.parse_args()
npzfile=args.f
mdl='_'.join(npzfile.split('/')[-1].split('_')[:2])
T=npzfile.split('/')[-1][:-4].split('_')[-1]
outname=f'cont_life/{mdl}_{T}'
dt=10e-6*10000 ### ns
Pcut=1.0

""""""
P=np.load(npzfile)['Pchn']
P=np.heaviside(P-Pcut,0.0)
nfr,nch=np.shape(P)[:2]

Ppad=np.pad(P,((1,1),(0,0),(0,0)),'constant',constant_values=-1)
dP=Ppad[1:]-Ppad[:-1]
Ps=np.where(dP[:-1]!=0,1,0)
Pe=np.where(dP[1:]!=0,1,0)

def logbin(exp1,exp2):
    bin=[]
    for i in range(exp1,exp2):
        for j in range(1,10):
            bin+=[j*10**i]
    return np.array(bin)
bin=logbin(-1,3)
Len0=[]
Len1=[]
# Len0log=[]
# Len1log=[]
for i in range(nch):
    for j in range(i):
        seglen=[]
        segval=[]
        idxs=np.where(Ps[:,i,j]==1)[0]
        idxe=np.where(Pe[:,i,j]==1)[0]
        for k,l in zip(idxs,idxe):
            seglen+=[l-k+1]
            segval+=[P[k,i,j]]
        seglen=np.array(seglen)
        segval=np.array(segval)
        len0=seglen[segval==0]
        len1=seglen[segval==1]
        len0=len0[len0>1]
        len1=len1[len1>1]
        Len0+=(len0*dt).tolist()
        Len1+=(len1*dt).tolist()
        # Len0log+=np.log10(np.array([x for x in len0 if x!=0])*dt).tolist()
        # Len1log+=np.log10(np.array([x for x in len1 if x!=0])*dt).tolist()
Len0=np.array(Len0)
Len1=np.array(Len1)
# Len0log=np.array(Len0log)
# Len1log=np.array(Len1log)
hist,edge=np.histogram(Len0,bins=bin)
# hist,edge=np.histogram(Len0log,bins=nhi)
T0=(edge[:-1]+edge[1:])/2
p_T0=hist/(edge[1:]-edge[:-1])/np.sum(hist)
T0_mean=np.exp(np.mean(np.log(Len0)))
T0_std=np.exp(np.std(np.log(Len0)))
Len0l=Len0[Len0<T0_mean]
T0_stdl=np.exp(np.sqrt(np.mean((np.log(Len0l)-np.log(T0_mean))**2)))
Len0r=Len0[Len0>T0_mean]
T0_stdr=np.exp(np.sqrt(np.mean((np.log(Len0r)-np.log(T0_mean))**2)))
hist,edge=np.histogram(Len1,bins=bin)
# hist,edge=np.histogram(Len1log,bins=nhi)
T1=(edge[:-1]+edge[1:])/2
p_T1=hist/(edge[1:]-edge[:-1])/np.sum(hist)
T1_mean=np.exp(np.mean(np.log(Len1)))
T1_std=np.exp(np.std(np.log(Len1)))
Len1l=Len1[Len1<T1_mean]
T1_stdl=np.exp(np.sqrt(np.mean((np.log(Len1l)-np.log(T1_mean))**2)))
Len1r=Len1[Len1>T1_mean]
T1_stdr=np.exp(np.sqrt(np.mean((np.log(Len1r)-np.log(T1_mean))**2)))

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.savez(f'{outname}.npz',T0=T0,p_T0=p_T0,T0_mean=T0_mean,T0_std=T0_std,T0_stdl=T0_stdl,T0_stdr=T0_stdr,
         T1=T1,p_T1=p_T1,T1_mean=T1_mean,T1_std=T1_std,T1_stdl=T1_stdl,T1_stdr=T1_stdr)
