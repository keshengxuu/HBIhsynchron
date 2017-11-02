# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 10:20:45 2015

@author: ksxuu
"""
from __future__ import print_function
import csv
import numpy as np
import matplotlib.pyplot as plt
tEnd=50000
nsim=50

def burst_detect(spikes,maxIntra=100,minInter=30):
    ISI=np.diff(spikes)
    num=len(ISI)
    
    ratio=2.5
    bflag=0
    Btime=[]
    Events=[]

    for i in range(num-1):
        if bflag==0:
            Btime.append(spikes[i])
            Events.append(1)
            if ((ISI[i]/ISI[i+1]>ratio and ISI[i+1]<maxIntra) or ISI[i+1]<minInter):
                bflag=1
        else:
            Events[-1]+=1
            if ((ISI[i+1]/ISI[i]>ratio and ISI[i+1]>minInter) or ISI[i+1]>maxIntra):
                bflag=0;
    CT=np.diff(Btime)                
    return Btime,CT,Events

