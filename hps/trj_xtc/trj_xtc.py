"""Import Modules"""
import argparse
import os
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
from trj import trj_rd
from xtc import xtc_wrt

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-f',type=str,help='The .trj file')
args=parser.parse_args()
trjfile=args.f
xtcfile=trjfile[:-4]+'.xtc'

""""""
r,tp,b=trj_rd(trjfile)
xtc_wrt(xtcfile,r,tp,b)