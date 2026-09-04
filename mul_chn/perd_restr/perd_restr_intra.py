"""Import Modules"""
import argparse
import numpy as np
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from xtc import xtc_rd,xtc_wrt
from perd_restr import str_restr

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-i',type=str,help='The inputting trajectory file')
parser.add_argument('-o',type=str,help='The outputting trajectory file')
args=parser.parse_args()
ixtcfile=args.i
oxtcfile=args.o
nat=163

""""""
r,tp,b=xtc_rd(ixtcfile,f'{ixtcfile[:-8]}.dat')
N=np.shape(r)[1]
nch=round(N/nat)
mol=[]
ref=[]
for i in range(nch):
    mol.append([i*nat,(i+1)*nat])
    ref+=[round(nat/2)]
r=str_restr(r,b,mol,ref)

""""""
xtc_wrt(oxtcfile,r,tp,b)





