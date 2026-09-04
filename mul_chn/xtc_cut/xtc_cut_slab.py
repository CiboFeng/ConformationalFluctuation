"""Import Modules"""
import argparse
import numpy as np
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from xtc import xtc_rd,xtc_wrt

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-i',type=str,help='The inputting trajectory file')
parser.add_argument('-o',type=str,help='The outputting trajectory file')
args=parser.parse_args()
ixtcfile=args.i
oxtcfile=args.o
sfr=5000
nfr=10000

""""""
r,tp,b=xtc_rd(ixtcfile,f'{ixtcfile[:-4]}.dat')
r=r[sfr:sfr+nfr]
xtc_wrt(oxtcfile,r,tp,b)
