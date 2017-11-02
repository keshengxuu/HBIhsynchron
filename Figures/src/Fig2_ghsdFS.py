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
#call the filter function
mle = fc.mle_filter(Table)
FR=fc.FR_filter(Table_FR)
FR_mask=fc.spikes_classfy(Table_FP)# firing pattern mask

#do the plot
Fig.Fig2_FS( mle, FR, FR_mask, tupTable2)


