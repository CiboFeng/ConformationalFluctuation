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
from xtc import xtc_rd,xtc_wrt

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-i',type=str,help='The inputing trajectory file')
parser.add_argument('-o',type=str,help='The outputing trajectory file')
args=parser.parse_args()
ixtcfile=args.i
oxtcfile=args.o
init=int(ixtcfile.split('/')[-1][:-4].split('_')[1])
rep=int(ixtcfile.split('/')[-1][:-4].split('_')[2])
nfrs=[1000,100,100,100,
      100,1000,100,100,
      100,100,100,1000,
      100,1000,100,1000]

""""""
r,tp,b=xtc_rd(ixtcfile,f'{ixtcfile[:-4]}.pdb')
r=r[nfrs[init]:]
xtc_wrt(oxtcfile,r,tp,b)

