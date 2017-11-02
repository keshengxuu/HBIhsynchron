# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:46:06 2016

@author: keshengXu
"""
import numpy as np

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
os.chdir(parentdir)


from lib import HBIh_model as HBM
from lib import Figures_plot as Fig


# do the simulation
Tstop = 20000
dt=0.05
time = np.arange(0,Tstop,dt)

gsd0=0.2846;gh0=0.1829 # for nonchaos
gsd1=0.203 ;gh1=0.335 # for chaos

Var_t = HBM.HBIh(time,gsd0,gh0)  # time series of Isolate nonchaotic neurons
Var_t1 = HBM.HBIh(time,gsd1,gh1) # time series of Isolate chaotic neurons

# do the plot
#Fig.Fig1_Schematic(time, Var_t1, Var_t)
Fig.Fig1_regluar(time, Var_t1, Var_t)





