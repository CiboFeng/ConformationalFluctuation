import argparse
import numpy as np
from scipy.spatial.distance import squareform,pdist
import os
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.figure import figaspect
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from xtc import xtc_rd

""""""
parser=argparse.ArgumentParser()
parser.add_argument('-d',type=str,help='The directory to trajectory files')
args=parser.parse_args()
xtcdir=args.d
outname='cont_prob_mom_1.5sgm'
mu=1
nmo=5
csgm=float(outname.split('_')[-1][:-3])
sgm=[5.04,6.56,5.68,5.58,5.48,6.02,5.92,4.50,6.08,6.18,6.18,6.36,6.18,6.36,5.56,5.18,5.62,6.78,6.46,5.86]
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)

""""""
r=np.zeros((0,nat,3))
xtcfiles=[f for f in os.listdir(xtcdir) if f.endswith('.xtc')]
for i in range(len(xtcfiles)):
    try:
        r=np.concatenate((r,xtc_rd(f'{xtcdir}/{xtcfiles[i]}',f'{xtcdir}/{xtcfiles[i][:-4]}.dat')[0]),axis=0)
    except Exception as err:
        print(err)

nfr=np.shape(r)[0]
Sgm=np.zeros(nat)
for i in range(nat):
    Sgm[i]=sgm[res.index(seq[i])]
Sgm=(Sgm[np.newaxis,:]+Sgm[:,np.newaxis])/2
P=np.zeros((nat,nat,nmo))
for i in range(nfr):
    P0=0.5*(1-np.tanh(mu*(squareform(pdist(r[i],'euclidean'))-csgm*Sgm)))
    for j in range(nmo):
        P[:,:,j]+=P0**(j+1)
P/=nfr

""""""
np.save(f'{outname}.npy',P)

