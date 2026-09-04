"""Import Modules"""
import argparse
import numpy as np
import os
import sys
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('D:\\Work\\Code\\Functions')
from tot_len import tot_len
from pyw import pyw
from xtc import xtc_rd
from perd_restr import str_restr

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-c',type=str,help='The .cfg file')
parser.add_argument('-t',type=str,help='The .trj file')
args=parser.parse_args()
cfgfile=args.c
trjfile=args.t
mdl='_'.join(trjfile.split('/')[-2].split('_')[-3:-1])
T=trjfile.split('/')[-2].split('_')[-1]
rawfilep=f'/hpc2hdd/home/cfeng593/proj/sim/rg_var/mul_chn/slab_{mdl}_{T}/fus'
outname=f'rerun_prs_{mdl}_{T}'
nthermo=1

"""Write to .cfg File"""
cfg=open(cfgfile,'w')
cfg.write(f'log\t\t{outname}.log\n')
cfg.write('units\t\treal\n')
cfg.write('boundary\tp p p\n')
cfg.write('timestep\t10.0\n')
cfg.write('dielectric\t\t80.0\n')
cfg.write('neighbor\t\t3.5 multi\n')
cfg.write('neigh_modify\tevery 10 delay 0\n')
cfg.write('special_bonds\tlj/coul 0.0 1.0 1.0\n\n')
cfg.write('atom_style\tfull\n')
cfg.write('bond_style\tharmonic\n')
cfg.write(f'pair_style\thybrid/overlay tanh/power 15.0 lj/lambda 0.1 0.0 35.0\n\n')
cfg.write(f'read_data\t{rawfilep}.dat\n')
cfg.write(f'include\t\t{rawfilep}.par\n\n')
cfg.write(f'thermo\t\t{nthermo}\n')
cfg.write('thermo_style\tcustom step pxx pyy pzz pxy pyz pxz\n\n')
cfg.write(f'rerun\t\t{trjfile} dump x y z wrapped no')
cfg.close()
