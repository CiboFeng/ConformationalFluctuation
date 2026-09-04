"""Import Modules"""
import argparse
import numpy as np
import os
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
from pyw import pyw
from tot_len import tot_len

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-f',type=str,help='The .log files')
args=parser.parse_args()
logfile=args.f
mdl='_'.join(logfile.split('/')[-2].split('_')[-2:])
idx=logfile.split('/')[-1].split('.')[0].split('_')[-1]
outname=f'mass_cent_dist/{mdl}_{idx}'
S=0.8

"""Read Log Data"""
d=[]
f=open(logfile,'r')
ls=f.readlines()
for i in range(len(ls)):
    l_s=ls[i].strip('\n').split()
    if len(l_s)==7:
        if l_s[0]=='Step' and l_s[6]=='v_dist':
            j=i+1
            while True:
                try:
                    l_s=ls[j].strip('\n').split()
                    d+=[float(l_s[6])]
                    j+=1
                except Exception:
                    break
            break
nfr=len(d)
nfr=round(nfr*S)
d=np.array(d[-nfr:])
z=np.abs((d-np.median(d))/(np.percentile(d,75)-np.percentile(d,25)))
d=d[z<10]

"""Output"""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
        except Exception:
            pass
f=open(f'{outname}.dat','w')
for i in range(nfr):
    f.write(f'{i} {d[i]}\n')
f.close()
