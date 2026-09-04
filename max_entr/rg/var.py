"""Import Modules"""
import sys
sys.path.append('D:\\Work\\Code\\Functions')
from pyw import pyw
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

"""Set Arguments"""
file0='../../all_atom/rg/rg.pyw'
file1='fus_wt_mom1_type_1.37sgm_0.6sgm/rg_itr20.pyw'
file2='fus_wt_mom1_type_1.55sgm_1.0sgm/rg_itr20.pyw'

"""Read Data and Calculate"""
q=np.array(pyw(file0,'Probability'))
rg0=q[0]
p_rg0=q[1]
q=np.array(pyw(file1,'Probability'))
rg1=q[0]
p_rg1=q[1]
q=np.array(pyw(file2,'Probability'))
rg2=q[0]
p_rg2=q[1]
q=np.array(pyw(file0,'Mean'))
rg0_mean=q[0]
q=np.array(pyw(file1,'Mean'))
rg1_mean=q[0]
q=np.array(pyw(file2,'Mean'))
rg2_mean=q[0]

""""""
var0=np.sqrt(np.sum(p_rg0*(rg0-rg0_mean)**2)*(rg0[1]-rg0[0]))
var1=np.sqrt(np.sum(p_rg1*(rg1-rg1_mean)**2)*(rg1[1]-rg1[0]))
var2=np.sqrt(np.sum(p_rg2*(rg2-rg2_mean)**2)*(rg2[1]-rg2[0]))
print(var0,var1,var2)