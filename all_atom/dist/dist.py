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
parser.add_argument('-d',type=str,help='The directory to trajectory files')
args=parser.parse_args()
xtcdir=args.d
outname=f'dist'
nat=163
nhi=50

""""""
r=np.zeros((0,nat,3))
xtcfiles=[f for f in os.listdir(xtcdir) if f.endswith('.xtc')]
for i in range(len(xtcfiles)):
    try:
        r=np.concatenate((r,xtc_rd(f'{xtcdir}/{xtcfiles[i]}',f'{xtcdir}/{xtcfiles[i][:-4]}.pdb')[0]),axis=0)
    except Exception as err:
        print(err)

nfr=np.shape(r)[0]
D=np.zeros((nfr,nat,nat))
for i in range(nfr):
    D[i]=squareform(pdist(r[i],'euclidean'))

d2=np.zeros((nat,nat,nhi))
p_d2=np.zeros((nat,nat,nhi))
d2_mean=np.mean(D,axis=0)
d2_std=np.std(D,axis=0)
d1=np.zeros((nat,nhi))
p_d1=np.zeros((nat,nhi))
d1_mean=np.zeros(nat)
d1_std=np.zeros(nat)
for i in range(nat):
    for j in range(i):
        hist,edge=np.histogram(D[:,i,j],bins=nhi)
        d2[i,j]=(edge[:-1]+edge[1:])/2
        d2[j,i]=d2[i,j]
        p_d2[i,j]=hist/np.sum(hist*(np.max(d2[i,j])-np.min(d2[i,j]))/nhi)
        p_d2[j,i]=p_d2[i,j]
    D_tmp=np.zeros((nfr,nat-i))
    for j in range(nat-i):
        D_tmp[:,j]=D[:,j,i+j]
    hist,edge=np.histogram(D_tmp,bins=nhi)
    d1[i]=(edge[:-1]+edge[1:])/2
    p_d1[i]=hist/np.sum(hist*(np.max(d1[i])-np.min(d1[i]))/nhi)
    d1_mean[i]=np.mean(D_tmp)
    d1_std[i]=np.std(D_tmp)

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.savez(f'{outname}.npz',d2=d2,p_d2=p_d2,d2_mean=d2_mean,d2_std=d2_std,d1=d1,p_d1=p_d1,d1_mean=d1_mean,d1_std=d1_std)