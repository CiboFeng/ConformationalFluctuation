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
from trj import trj_wrt

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-i',type=str,help='The inputting xtc file')
parser.add_argument('-o',type=str,help='The outputting trj file')
args=parser.parse_args()
xtcfile=args.i
trjfile=args.o
step=int(1e4)

""""""
r,tp,b=xtc_rd(xtcfile,f'{xtcfile[:-8]}.dat')
b=b-np.mean(b,axis=-1,keepdims=True)
trj_wrt(trjfile,r,tp,b,step=step)