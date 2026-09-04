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
outname='rg'
nat=163
nhi=200

""""""
r=np.zeros((0,nat,3))
xtcfiles=[f for f in os.listdir(xtcdir) if f.endswith('.xtc')]
for i in range(len(xtcfiles)):
    try:
        r=np.concatenate((r,xtc_rd(f'{xtcdir}/{xtcfiles[i]}',f'{xtcdir}/{xtcfiles[i][:-4]}.pdb')[0]),axis=0)
    except Exception as err:
        print(err)
rc=np.mean(r,axis=1,keepdims=True)
rg=np.mean(np.linalg.norm(r-rc,axis=-1),axis=-1)
rg_mean=np.mean(rg)
rg_std=np.std(rg)
hist,edge=np.histogram(rg,bins=nhi)
rg=(edge[:-1]+edge[1:])/2
p_rg=hist/np.sum(hist*(np.max(rg)-np.min(rg))/nhi)

""""""
f=open(f'{outname}.pyw','w')
f.write(f'Radius of Gyration, Its Probability Density:\n')
for i in range(nhi):
    f.write(f'{rg[i]} {p_rg[i]}\n')
f.write('\n')
f.write('Mean of Radius of Gyration:\n')
f.write(f'{rg_mean}\n\n')
f.write('Standard Deviation of Radius of Gyration:\n')
f.write(f'{rg_std}')
f.close()

