"""Import Modules"""
import argparse
import numpy as np
import os
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from xtc import xtc_rd,xtc_wrt
from perd_restr import str_restr_dyn

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-f',type=str,help='The trajectory file')
args=parser.parse_args()
xtcfile=args.f
mdl='_'.join(xtcfile.split('/')[-2].split('_')[-3:-1])
T=xtcfile.split('/')[-2].split('_')[-1]
outname=f'end_dist_chk_slab/{mdl}_{T}'
nat=163
nhi=100

""""""
r,tp,b=xtc_rd(xtcfile,f'{xtcfile[:-14]}.dat')
nfr,N,_=np.shape(r)
nch=round(N/nat)
d=np.zeros((nfr,nch,4))
r=r.reshape(nfr,nch,nat,3)
dr=r[:,:,0]-r[:,:,-1]
d[:,:,:3]=np.abs(dr)
d[:,:,3]=np.linalg.norm(dr,axis=-1)
D=np.zeros((4,nhi))
p_D=np.zeros((4,nhi))
for i in range(4):
    hist,edge=np.histogram(d[:,:,i],bins=nhi)
    D[i]=(edge[:-1]+edge[1:])/2
    p_D[i]=hist/np.sum(hist*(np.max(D[i])-np.min(D[i]))/nhi)

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
pyw=open(f'{outname}.pyw','w')
pyw.write('End-end Distance, Its Probability Density of x Direction, Similar for y Direction, z Direction, and All:\n')
for i in range(nhi):
    for j in range(4):
        pyw.write(f'{D[j,i]} {p_D[j,i]} ')
    pyw.write('\n')
pyw.close()