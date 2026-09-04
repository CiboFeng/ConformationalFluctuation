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
mdl='_'.join(xtcfile.split('/')[-2].split('_')[-3:-1])
T=xtcfile.split('/')[-2].split('_')[-1]
outname=f'jump/{mdl}_{T}'
res=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
m=[71.08,156.20,114.10,115.10,103.10,128.10,129.10,57.05,137.10,113.20,113.20,128.20,131.20,147.20,97.12,87.08,101.10,186.20,163.20,99.07]
seq='MASNDYTQQA TQSYGAYPTQ PGQGYSQQSS QPYGQQSYSG YSQSTDTSGY GQSSYSSYGQ SQNTGYGTQS TPQGYGSTGG YGSSQSSQSS YGQQSSYPGY ' \
    'GQQPAPSSTS GSYGSSSQSS SYGQPQSGSY SQQPSYGGQQ QSYGQQQSYN PPQGYGQQNQ YNS'
seq=seq.replace(' ','')
nat=len(seq)
seq=[res.index(seq[i]) for i in range(nat)]
nhi=200

""""""
r=xtc_rd(xtcfile,f'{xtcfile[:-12]}.dat')[0]
nfr,N,_=np.shape(r)
nch=round(N/nat)
r=r.reshape(nfr,nch,nat,3)
m=np.array([m[seq[i]] for i in range(nat)])
M=np.tile(m,(nch,1))
rc=np.sum(M[np.newaxis,:,:,np.newaxis]*r,axis=2)/np.sum(m)

d=np.linalg.norm(r[1:]-r[:-1],axis=-1)
D=np.zeros((nat,nhi))
p_D=np.zeros((nat,nhi))
for i in range(nat):
    hist,edge=np.histogram(d[:,:,i],bins=nhi)
    D[i]=(edge[:-1]+edge[1:])/2
    p_D[i]=hist/np.sum(hist*(np.max(D[i])-np.min(D[i]))/nhi)
dc=np.linalg.norm(rc[1:]-rc[:-1],axis=-1)
hist,edge=np.histogram(dc,bins=nhi)
Dc=(edge[:-1]+edge[1:])/2
p_Dc=hist/np.sum(hist*(np.max(Dc)-np.min(Dc))/nhi)

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
pyw=open(f'{outname}.pyw','w')
pyw.write('Displacement, It Probability Density, Similar for Each Residue:\n')
for i in range(nhi):
    for j in range(nat):
        pyw.write(f'{D[j,i]} {p_D[j,i]} ')
    pyw.write(f'\n')
pyw.write('\n')
pyw.write('Displacement of Center of Mass, It Probability Density:\n')
for i in range(nhi):
    pyw.write(f'{Dc[i]} {p_Dc[i]}\n')
pyw.close()