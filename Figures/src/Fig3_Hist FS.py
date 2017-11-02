# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:46:06 2016

@author: keshengxuu@gmail.com
"""
import numpy as np
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
os.chdir(parentdir)


from lib import loaddata as ld
from lib import filter_function as fc
from lib import Figures_plot as Fig
execfile('lib/Defaultsets.py')

#load the data
tupTable1, tupTable2  = ld.Fig2_load()
Table, Table_FP, Table_FR = tupTable1
sp_chao,sp_nonchao=ld.load_data_spike(rseed0,ranges0)
#for i in range(5):
Burst_chaos, Burst_nonch = fc.frac_burst(sp_chao,sp_nonchao, Nindex)


chaosTableFR=Table_FR[(Table_FR[:,0]>ghmin1/1000)*(Table_FR[:,0]<ghmax1/1000)*
                    (Table_FR[:,1]>gsdmin1/1000)*(Table_FR[:,1]<gsdmax1/1000)]
nonchaosTableFR=Table_FR[(Table_FR[:,0]>ghmin2/1000)*(Table_FR[:,0]<ghmax2/1000)*
                    (Table_FR[:,1]>gsdmin2/1000)*(Table_FR[:,1]<gsdmax2/1000)]
chaosTableFR2=Table_FR[(Table_FR[:,0]>ghmin1b/1000)*(Table_FR[:,0]<ghmax1b/1000)*
                    (Table_FR[:,1]>gsdmin1b/1000)*(Table_FR[:,1]<gsdmax1b/1000)]
nonchaosTableFR2=Table_FR[(Table_FR[:,0]>ghmin2b/1000)*(Table_FR[:,0]<ghmax2b/1000)*
                    (Table_FR[:,1]>gsdmin2b/1000)*(Table_FR[:,1]<gsdmax2b/1000)]

tupTable = (chaosTableFR, nonchaosTableFR, chaosTableFR2, nonchaosTableFR2 )
#do the plot
Fig.Fig3Histgram_FS(tupTable, Burst_chaos, Burst_nonch)


