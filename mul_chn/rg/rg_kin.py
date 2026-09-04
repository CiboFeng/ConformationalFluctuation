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
mdl='_'.join(xtcfile.split('/')[-2].split('_')[-2:])
rep=xtcfile.split('/')[-1][:-14].split('_')[1]
outname=f'rg_kin/{mdl}_{rep}'
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]

""""""
m=np.array([m[seq[i]] for i in range(nat)])
r=xtc_rd(xtcfile,f'{xtcfile[:-14]}.dat')[0]
nfr,N,_=np.shape(r)
nch=round(N/nat)
r=r.reshape(nfr,nch,nat,3)
rg=np.zeros((nfr,nch))
rc=np.sum(m[np.newaxis,np.newaxis,:,np.newaxis]*r,axis=2,keepdims=True)/np.sum(m)
rg=np.sum(m[np.newaxis,np.newaxis,:]*np.linalg.norm(r-rc,axis=-1),axis=-1)/np.sum(m)

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
pyw=open(f'{outname}.pyw','w')
pyw.write(f'Radius of Gyration of the 0-{nch-1}-th Chain:\n')
for i in range(nfr):
    for j in range(nch):
        pyw.write(f'{rg[i,j]} ')
    pyw.write('\n')
pyw.close()
