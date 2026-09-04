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
init=int(xtcfile.split('/')[-1][:-4].split('_')[1])
rep=int(xtcfile.split('/')[-1][:-4].split('_')[2])
outname=f'rmsd/{init}_{rep}'

""""""
r=xtc_rd(xtcfile,f'{xtcfile[:-4]}.pdb')[0]
nfr=np.shape(r)[0]
D0=squareform(pdist(r[0],'euclidean'))
d=np.zeros(nfr)
for i in range(nfr):
    d[i]=(np.mean((squareform(pdist(r[i],'euclidean'))-D0)**2))**0.5

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
f=open(f'{outname}.pyw','w')
f.write(f'Index of Frames, Root Mean Square Deviation:\n')
for i in range(nfr):
    f.write(f'{i} {d[i]}\n')
f.close()

