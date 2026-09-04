"""Import Modules"""
import argparse
import numpy as np
import os
import random
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/opt/mypylib')
sys.path.append('/hpc2hdd/home/chuo657/cbfengphy/opt/mypylib')
from pyw import pyw
from tot_len import tot_len

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-s',type=str,help='The sequence file')
parser.add_argument('-c',type=str,help='The written coefficient file')
parser.add_argument('-d',type=str,help='The .dat file prefix')
parser.add_argument('-nm',type=int,help='The number of moment')
parser.add_argument('-np',type=int,help='The number of parallel trajectories')
args=parser.parse_args()
seqfile=args.s
coefffile=args.c
datpre=args.d
nmo=args.nm
n=args.np
ntp=20
R0=3.82
dspace=10.0
chrg=[0,1,0,-1,0,0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0]

"""Read Sequence Data and Determine the System Size"""
seq=pyw(seqfile,'Sequence','int')[0]
N=len(seq)

"""Write to Initial Coefficient File"""
pyw=open(coefffile,'w')
pyw.write(f'The 1-{nmo}th Coefficient of Data-driven Term:\n')
for i in range(ntp):
    for j in range(ntp):
        for k in range(nmo):
            pyw.write(f'{0.0}\t')
        pyw.write('\n')
pyw.close()

"""Write to .dat Files Cyclically"""
for i in range(0,n):
    #"""Generate an Initial Conformation"""
    x=[0]
    y=[0]
    z=[0]
    for j in range(1,N):
        judge=True
        while judge:
            theta=random.uniform(0,np.pi)
            phi=random.uniform(0,2*np.pi)
            x1=x[j-1]+R0*np.sin(theta)*np.cos(phi)
            y1=y[j-1]+R0*np.sin(theta)*np.sin(phi)
            z1=z[j-1]+R0*np.cos(theta)
            dist=[]
            for k in range(len(x)):
                dist+=[((x1-x[k])**2+(y1-y[k])**2+(z1-z[k])**2)**0.5]
            dmin=min(dist)
            if dmin>R0:
                judge=False
        x+=[x1]
        y+=[y1]
        z+=[z1]

    #"""Calculate the Spacial Range of the Box"""
    xlo=min(x)-dspace*(max(x)-min(x))
    xhi=max(x)+dspace*(max(x)-min(x))
    ylo=min(x)-dspace*(max(y)-min(y))
    yhi=max(x)+dspace*(max(y)-min(y))
    zlo=min(x)-dspace*(max(z)-min(z))
    zhi=max(x)+dspace*(max(z)-min(z))

    #"""Write to .dat File"""
    dat=open(f'{datpre}_{i}.dat','w')
    dat.write('LAMMPS Data\n\n')
    dat.write(f'{N}\tatoms\n')
    dat.write(f'{N-1}\tbonds\n\n')
    dat.write(f'{ntp}\tatom types\n')
    dat.write(f'1\tbond types\n\n')
    dat.write(f'{tot_len(xlo,16)}\t{tot_len(xhi,16)}\txlo\txhi\n')
    dat.write(f'{tot_len(ylo,16)}\t{tot_len(yhi,16)}\tylo\tyhi\n')
    dat.write(f'{tot_len(zlo,16)}\t{tot_len(zhi,16)}\tzlo\tzhi\n\n')
    dat.write('Atoms\n\n')
    for j in range(N):
        dat.write(f'{j+1}\t1\t{seq[j]+1}\t{tot_len(chrg[seq[j]],16)}\t{tot_len(x[j],16)}\t{tot_len(y[j],16)}\t{tot_len(z[j],16)}\n')
    dat.write('\n')
    dat.write('Bonds\n\n')
    for j in range(N-1):
        dat.write(f'{j+1}\t1\t{j+1}\t{j+2}\n')
    dat.write('\n')
    dat.close()