#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May 20 19:15:27 2017

@author: ksxuu
"""
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
os.chdir(parentdir)

from lib import loaddata as ld
from lib import Figures_plot as Fig
from lib import filter_function as fc


# load the data 
chaos, nonchaos, noIh = ld.Fig5_load()

# do the plot
Fig.Fig5_plot(chaos, nonchaos, noIh)


