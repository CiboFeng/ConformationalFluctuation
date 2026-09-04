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
from xtc import xtc
from trj import trj

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-d',type=str,help='The directory to trajectory files')
args=parser.parse_args()
xtcdir=args.d
itr=xtcdir.split('_')[-1]
outname=f'rg_{itr}'
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
nhist=50

""""""
seq=[res.index(seq[i]) for i in range(nat)]
m=np.array([m[seq[i]] for i in range(nat)])

r=np.zeros((0,nat,3))
xtcfiles=[f for f in os.listdir(xtcdir) if f.endswith('.xtc')]
for i in range(len(xtcfiles)):
    try:
        r=np.concatenate((r,xtc(f'{xtcdir}/{xtcfiles[i]}',f'{xtcdir}/{xtcfiles[i][:-4]}.dat')[-1]),axis=0)
    except:
        pass
rc=np.sum(m[np.newaxis,:,np.newaxis]*r,axis=1,keepdims=True)/np.sum(m)
rg=np.sum(m[np.newaxis,:]*np.linalg.norm(r-rc,axis=-1),axis=-1)/np.sum(m)
rg_mean=np.mean(rg)
hist,edge=np.histogram(rg,bins=nhist)
rg=(edge[:-1]+edge[1:])/2
p_rg=hist/np.sum(hist*(np.max(rg)-np.min(rg))/nhist)

""""""
f=open(f'{outname}.pyw','w')
f.write(f'Radius of Gyration, Its Probability Density:\n')
for i in range(nhist):
    f.write(f'{rg[i]} {p_rg[i]}\n')
f.write('\n')
f.write('Mean of Radius of Gyration:\n')
f.write(f'{rg_mean}')
f.close()