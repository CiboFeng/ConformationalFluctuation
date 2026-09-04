"""Import Modules"""
import argparse
import numpy as np
import networkx as nx
import os
import sys
sys.path.append('D:\\Work\\Code\\Functions')
sys.path.append('/hpc2hdd/home/cfeng593/opt/mypylib')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from pyw import pyw

""""""
parser=argparse.ArgumentParser()
parser.add_argument('-f',type=str,help='The .npz file')
args=parser.parse_args()
npzfile=args.f
mdl='_'.join(npzfile.split('/')[-1].split('_')[:2])
T=npzfile.split('/')[-1][:-4].split('_')[-1]
outname=f'netw/{mdl}_{T}'
nhi=20
eul=0.5772156649015328606065120900824024310421
Pcut=1.0

""""""
P=np.load(npzfile)['Pchn']
P=np.heaviside(P-Pcut,0.0)
nfr,nch=np.shape(P)[:2]
qnode=np.zeros((nfr,nch,2))
qnetw=np.zeros((nfr,9))
for i in range(nfr):
    G=nx.from_numpy_array(P[i])

    nedge=G.number_of_edges()
    deg=list(dict(G.degree()).values())
    deg_mean=np.mean(deg)

    betweenness=nx.betweenness_centrality(G,normalized=True)
    bc=np.array(list(betweenness.values()))
    bc_mean=np.mean(bc)
    bc_std=np.std(bc)
    bc_skwn=np.mean((bc-bc_mean)**3)/bc_std**3
    # threshold=np.percentile(bc,95)
    # hub=[node for node,bc in betweenness.items() if bc > threshold]

    clst=nx.average_clustering(G)
    clst0=2*nedge/(nch*(nch-1))
    if nx.is_connected(G):
        lpath=nx.average_shortest_path_length(G)
        lpath0=(np.log(nch)-eul)/np.log(2*nedge/nch)+0.5
    else:
        largest_cc=max(nx.connected_components(G),key=len)
        G_lcc=G.subgraph(largest_cc)
        lpath=nx.average_shortest_path_length(G_lcc)
        nnode_lcc=G_lcc.number_of_nodes()
        nedge_lcc=G_lcc.number_of_edges()
        lpath0=(np.log(nnode_lcc)-eul)/np.log(2*nedge_lcc/nnode_lcc)+0.5
    clst_rel=clst/clst0 
    lpath_rel=lpath/lpath0
    sml=clst_rel/lpath_rel ### Small-world coefficient (Humphries & Gurney, 2008).

    qnode[i,:,0]=deg
    qnode[i,:,1]=bc
    qnetw[i,0]=deg_mean
    qnetw[i,1]=bc_mean
    qnetw[i,2]=bc_std
    qnetw[i,3]=bc_skwn
    qnetw[i,4]=clst
    qnetw[i,5]=lpath
    qnetw[i,6]=clst_rel
    qnetw[i,7]=lpath_rel
    qnetw[i,8]=sml

x=np.zeros((nhi,11))
p=np.zeros((nhi,11))
x_mean=np.zeros(11)
x_std=np.zeros(11)
for i in range(2):
    hist,edge=np.histogram(qnode[:,:,i],bins=nhi)
    x[:,i]=(edge[:-1]+edge[1:])/2
    p[:,i]=hist/np.sum(hist*(np.max(x[:,i])-np.min(x[:,i]))/nhi)
    x_mean[i]=np.mean(qnode[:,:,i])
    x_std[i]=np.std(qnode[:,:,i])
for i in range(9):
    hist,edge=np.histogram(qnetw[:,i],bins=nhi)
    x[:,i+2]=(edge[:-1]+edge[1:])/2
    p[:,i+2]=hist/np.sum(hist*(np.max(x[:,i+2])-np.min(x[:,i+2]))/nhi)
    x_mean[i+2]=np.mean(qnetw[:,i])
    x_std[i+2]=np.std(qnetw[:,i])

""""""
if '/' in outname:
    outdir='/'.join(outname.split('/')[:-1])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
np.savez(f'{outname}.npz',x=x,p=p,x_mean=x_mean,x_std=x_std)