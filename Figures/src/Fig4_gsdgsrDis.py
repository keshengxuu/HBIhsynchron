# -*- coding: utf-8 -*-
"""
Firing rate statistics
FIXED DISTRIBUTION OF gsd, gh

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
Table3, Table4 = ld.FIg4_load()
#set  default  maxima and minima value of firing rate
minFR=3
maxFR=4.5
sim=50
# choose chao, nonchaotic and noIh table,respectively.
chaosTable=Table3[(Table3[:,2]>minFR)*(Table3[:,2]<maxFR)*(Table3[:,-1]>0.001)]
nonchaosTable=Table3[(Table3[:,2]>minFR)*(Table3[:,2]<maxFR)*(Table3[:,-1]<0.00008)]
nonIhTable=Table4[(Table4[:,2]>minFR)*(Table4[:,2]<maxFR)*(Table4[:,-1]<0.000001)*(Table4[:,-1]>-0.001)]
nonchaosTable2=[]
nonIhTable2=[]
for freq in chaosTable[:,2]:
    if freq in nonchaosTable[:,2]:
        indices=np.where(nonchaosTable[:,2]==freq)[0]
        nonchaosTable2.append(nonchaosTable[np.random.choice(np.where(nonchaosTable[:,2]==freq)[0],1)[0]])
    else:
        nonchaosTable2.append(nonchaosTable[np.abs(nonchaosTable[:,2]-freq).argmin()])

    if freq in nonIhTable[:,2]:
        indices=np.where(nonIhTable[:,2]==freq)[0]
        nonIhTable2.append(nonIhTable[np.random.choice(np.where(nonIhTable[:,2]==freq)[0],1)[0]])
    else:
        nonIhTable2.append(nonIhTable[np.abs(nonIhTable[:,2]-freq).argmin()])

nonchaosTable2=np.array(nonchaosTable2)
nonIhTable2=np.array(nonIhTable2)

chaosindex=np.random.choice(len(chaosTable),size=sim,replace=False)
chaosvalues=chaosTable[chaosindex,0:4]
nonchaosindex=np.random.choice(len(nonchaosTable2),size=sim,replace=False)
nonchaosvalues=nonchaosTable2[nonchaosindex,0:4]
nonIhindex=np.random.choice(len(nonIhTable2),size=sim,replace=False)
nonIhvalues=nonIhTable2[nonIhindex,0:4]

# choosing the same distribution for chaos ,nonchaos and nonchaos without ih current
lyap=fc.mle_filter(Table3)
lyap2=fc.mle_filter(Table4)

Freq2=lyap[0]
Lyap2=lyap[1]
#mask chaos and nonchaos
chaosMask=0.5 * np.ones_like(Freq2) + 0.5* ((Freq2>minFR)*(Freq2<maxFR)*(Lyap2>0.001))
nonchaosMask=0.5 * np.ones_like(Freq2) + 0.5* ((Freq2>minFR)*(Freq2<maxFR)*(Lyap2<0.000))

# mask without ih current
Freqnoih=lyap2[0]
Lyapnoih=lyap2[1]
noihMask=0.5 * np.ones_like(Freqnoih) + 0.5* ((Freqnoih>minFR)*(Freqnoih<maxFR)*(Lyapnoih<0.000))

# maks the desired regions
auxchao=fc.mask_function(lyap,chaosMask)
auxnochaos=fc.mask_function(lyap,nonchaosMask)
auxnoih=fc.mask_function(lyap2,noihMask)

tupTable = (chaosTable, nonchaosTable2, nonIhTable2 )
# do the plot
Fig.Fig4_Dis(lyap, lyap2, auxchao, auxnochaos, auxnoih, tupTable)

