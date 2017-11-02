# -*- coding: utf-8 -*-
"""
Created on 20.April 17 22:38:58 2016
keshengxuu@gmail.com
@author: ksxuu

load_data_SBE's function :
INPUT
rseed : different seeds for simulation in the neural networks
ranges:  four different ranges  which i made

OUTPUT:
sbe_chaos: spikes trains of chaos from different ranges
sbe_nonchaos: soikes trains of nonchaos from different ranges


load_data_spikes()'s function:

table is as following:
Ggj,  np.mean(frInd),  np.std(frInd), CVspikes,  np.mean(phaseCoh),  CVenv
 0      1               2               3               4               5
 
chaotic :output the table of chaos
nonchaotic output the table if  nonchaos

"""
from __future__ import print_function
import csv
import numpy as np
import matplotlib.pyplot as plt



def load_data_spike(rseed,ranges):
    datapath = 'data/FigX_data/spiketrain/'
    def read_csv_file(filename):
        """Reads a CSV file and return it as a list of rows."""
        data = []
        for row in csv.reader(open(filename)):
    # transfer the string type to the float type
            row=list((float(x) for x in row))
            data.append(row)
        return data
    
    Ggjvals=np.logspace(-6,0,num=23)
    Ggjvals=np.insert(Ggjvals,0,0)
    
    path_chaos=[]
    path_nonchaos=[]
    sbe_chaos=[]
    sbe_nonchaos=[]
    for i,rs in zip(ranges,rseed):
        pat0=datapath+'chaos50nodes-seed%g/Spikes_s%g_G'%(rs,rs)
        pat1=datapath+'nonchaos50nodes-seed%g/Spikes_s%g_G'%(rs,rs)
        path_chaos.append(pat0)
        path_nonchaos.append(pat1)
        
    for pa0,pa1 in zip(path_chaos,path_nonchaos):
        sp0=[read_csv_file(pa0+'%.1e.txt'%s) for s in Ggjvals]
        sp1=[read_csv_file(pa1+'%.1e.txt'%s) for s in Ggjvals]
        sbe_chaos.append(sp0)
        sbe_nonchaos.append(sp1)
    return sbe_chaos,sbe_nonchaos
    


"""
Ggj,  np.mean(frInd),  np.std(frInd), CVspikes,  np.mean(phaseCoh),  CVenv
 0      1               2               3               4               5
 
"""
#def load_data_table():
#    #load the data from datafile
#    chaotic=np.array([np.loadtxt('data/chaotic'+'/Data_chaoticrange%g.txt'%s,delimiter=',') for s in range(1,5)])
#    nonchaotic=np.array([np.loadtxt('data/nonchaotic'+'/Data_nonchaoticrange%g.txt'%s,delimiter=',') for s in range(1,5)])
#    return chaotic,nonchaotic
    
def load_data_table(ranges):
    #load the data from datafile
    chaotic=np.array([np.loadtxt('data/chaotic'+'/Data_chaoticrange%g.txt'%s,delimiter=',') for s in ranges])
    nonchaotic=np.array([np.loadtxt('data/nonchaotic'+'/Data_nonchaoticrange%g.txt'%s,delimiter=',') for s in ranges])
    return chaotic,nonchaotic


def Fig2_load(Arguments=None):
    #load mle ,firing pattern and firing rate data of isolate neruons
    Table=np.loadtxt('data/Fig2_data/mlefinalTable_gsdgh36.txt',skiprows=1)
    Table_FP = np.loadtxt("data/Fig2_data/finalTable_ghgsd36.txt",skiprows=1)
    Table_FR=np.loadtxt("data/Fig2_data/finalTablegsdgh.txt",skiprows=1)
    #lood the synronization data 
    chaor1= np.loadtxt('data/Fig2_data/ChaosRange1-50nodes_table.txt',delimiter=',')
    nonchaor1= np.loadtxt('data/Fig2_data/NonchaosRange1-50nodes_table.txt',delimiter=',')
    chaor2= np.loadtxt('data/Fig2_data/ChaosRange2-50nodes_table.txt',delimiter=',')
    nonchaor2= np.loadtxt('data/Fig2_data/NonchaosRange2-50nodes_table.txt',delimiter=',')
    
    tupTable1 = (Table, Table_FP, Table_FR)
    tupTable2 = (chaor1,nonchaor1, chaor2, nonchaor2)
    return  tupTable1, tupTable2

def FIg4_load(Arguments=None):
    #load the datatxt
    Table3=np.loadtxt("data/Fig4_data/mleFR_gsdgsrzoomT36gh04.txt",skiprows=1,delimiter=',')
    Table4=np.loadtxt("data/Fig4_data/mleFR_gsdgsrzoomT36gh00_noih.txt",skiprows=1,delimiter=',')
    return Table3, Table4

def Fig5_load(Arguments=None):
    datapath = 'data/Fig5_data/'
    FILENA=['FR30to45','FR70to95']
    chaos=[np.loadtxt(datapath+'%schaos-50nodes_table.txt'%s,delimiter=',')for s in FILENA]
    nonchaos=[np.loadtxt(datapath+'%snonchaos-50nodes_table.txt'%s,delimiter=',') for s in  FILENA]
    noIh=[np.loadtxt(datapath+'%snoIh-50nodes_table.txt'%s,delimiter=',') for s in  FILENA]
    return chaos, nonchaos, noIh