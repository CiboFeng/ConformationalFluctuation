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
outname=f'comm_cont_prob/{idx}'
mu=1.0
csgm=1.5
sgm=[5.04,6.56,5.68,5.58,5.48,6.02,5.92,4.50,6.08,6.18,6.18,6.36,6.18,6.36,5.56,5.18,5.62,6.78,6.46,5.86]
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

"""Calculate"""
r=np.zeros((0,nat,3))
xtcfiles=[f for f in os.listdir(xtcdir) if f.endswith('.xtc')]
for i in range(len(xtcfiles)):
    try:
        r=np.concatenate((r,xtc_rd(f'{xtcdir}/{xtcfiles[i]}',f'{xtcdir}/{xtcfiles[i][:-4]}.dat')[0]),axis=0)
    except Exception as err:
        print(err)
nfr=np.shape(r)[0]

seq=[res.index(seq[i]) for i in range(nat)]
Sgm=np.zeros(nat)
for i in range(nat):
    Sgm[i]=sgm[seq[i]]
Sgm=(Sgm[np.newaxis,:]+Sgm[:,np.newaxis])/2
P=np.zeros((nat,nat))
for i in range(nfr):
    D=squareform(pdist(r[i],'euclidean'))
    P+=(0.5*(1-np.tanh(mu*(D-csgm*Sgm))))
P/=nfr

"""Output Data"""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
pyw=open(f'{outname}.pyw','w')
pyw.write('Contact Probability:\n')
for i in range(nat):
    for j in range(nat):
        pyw.write(f'{P[i,j]}\n')
pyw.close()