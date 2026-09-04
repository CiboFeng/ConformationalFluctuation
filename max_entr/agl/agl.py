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
idx='_'.join(xtcdir.split('_')[-3:])
outname=f'agl/{idx}'
nat=163
nhi=100

""""""
r=np.zeros((0,nat,3))
xtcfiles=[f for f in os.listdir(xtcdir) if f.endswith('.xtc')]
for i in range(len(xtcfiles)):
    try:
        r=np.concatenate((r,xtc_rd(f'{xtcdir}/{xtcfiles[i]}',f'{xtcdir}/{xtcfiles[i][:-4]}.dat')[0]),axis=0)
    except Exception as err:
        print(err)
a=r[:,:-2]-r[:,1:-1]
b=r[:,2:]-r[:,1:-1]
tht=np.arccos(np.sum(a*b,axis=-1)/np.linalg.norm(a,axis=-1)/np.linalg.norm(b,axis=-1))
Tht_mean=np.mean(tht,axis=0)
Tht_std=np.std(tht,axis=0)
Tht=np.zeros((nat-2,nhi))
p_Tht=np.zeros((nat-2,nhi))
bin=np.linspace(0,np.pi,nhi+1) ### linspace contains the end.
for i in range(nat-2):
    hist,edge=np.histogram(tht[:,i],bins=bin)
    Tht[i]=(edge[:-1]+edge[1:])/2
    p_Tht[i]=hist/np.sum(hist*np.pi/nhi)
Tht_mean_all=np.mean(tht)
Tht_std_all=np.std(tht)
hist,edge=np.histogram(tht,bins=bin)
Tht_all=(edge[:-1]+edge[1:])/2
p_Tht_all=hist/np.sum(hist*np.pi/nhi)

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
f=open(f'{outname}.pyw','w')
f.write(f'Each-residue Angle, Its Probability Density:\n')
for i in range(nat-2):
    for j in range(nhi):
        f.write(f'{Tht[i,j]} {p_Tht[i,j]}\n')
f.write('\n')
f.write('Mean of Each-residue Angle:\n')
for i in range(nat-2):
    f.write(f'{Tht_mean[i]}\n')
f.write('\n')
f.write('Standard Deviation of Each-residue Angle:\n')
for i in range(nat-2):
    f.write(f'{Tht_std[i]}\n')
f.write('\n')
f.write(f'All Angle, Its Probability Density:\n')
for i in range(nhi):
    f.write(f'{Tht_all[i]} {p_Tht_all[i]}\n')
f.write('\n')
f.write('Mean of All Angle:\n')
f.write(f'{Tht_mean_all}\n\n')
f.write('Standard Deviation of All Angle:\n')
f.write(f'{Tht_std_all}')
f.close()