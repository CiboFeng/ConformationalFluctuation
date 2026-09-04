"""Import Modules"""
import argparse
import datetime
import numpy as np
from scipy.spatial.distance import squareform,pdist
import MDAnalysis as mda
import sys
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
from xtc import xtc_rd,xtc_wrt

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument('-t',type=str,help='The topology file')
parser.add_argument('-i',type=str,help='The input .xtc file')
parser.add_argument('-o',type=str,help='The output .xtc file')
args=parser.parse_args()
topfile=args.t
ixtcfile=args.i
oxtcfile=args.o
pdbfile=oxtcfile[:-4]+'.pdb'

"""Read Data"""
print('Begin to read .xtc file:',datetime.datetime.now())
R,tp,b=xtc_rd(ixtcfile,topfile)
print('Finish reading .xtc file:',datetime.datetime.now())
u=mda.Universe(topfile)
atoms=[]
for atom in u.atoms:
    atoms.append({
        'atom_id': atom.id,          
        'atom_name': atom.name,   
        'atom_type': atom.type,
        'resid': atom.resid,
        'resname': atom.resname,
        'segment': atom.segment.segid
    })

"""Reduce"""
nfr,nat=np.shape(R)[:2]
r=np.zeros((nfr,0,3))
for i in range(nat):
    if atoms[i]['atom_name']=='CA':
        r=np.concatenate((r,R[:,i:i+1]),axis=1)
r=r[::10]
nfr,nat=np.shape(r)[:2]

"""Create New Topology and .pdb File"""
u=mda.Universe.empty(n_atoms=nat,n_residues=nat,trajectory=True,atom_resindex=list(range(nat)),residue_segindex=[0]*nat)
u.add_TopologyAttr('name',['CA']*nat)
u.add_TopologyAttr('type',['CA']*nat)
u.add_TopologyAttr('resnames',[tp[i] for i in range(nat)])
u.add_TopologyAttr('resids',list(range(1,nat+1)))
zero_coords=np.zeros((nat,3))
u.load_new(zero_coords,format='memory')
u.atoms.write(pdbfile)

"""Write to .xtc File"""
xtc_wrt(oxtcfile,r,tp,b)
