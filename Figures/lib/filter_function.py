# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 16:55:42 2015

@author: keshengxuu
"""
import numpy as np
import contextlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import rc, cm
import matplotlib as mlp
from numpy.ma import masked_array
from matplotlib.colors import Normalize
import StringIO

def spikes_classfy(Table):
    ghN=len(np.unique(Table[:,0]))
    gsdN=len(np.unique(Table[:,1]))
    a,b,c,d,e=[np.zeros([ghN, gsdN]) for i in range(5)]
    k=0
    #for i in range(0,91):
    #    for j in range(0,91):
    for i in range(0,ghN):
        for j in range(0,gsdN):
            c[i,j] = Table[k,6]
            if (Table[k,6]==4)and(Table[k,7]>100): #reconvierte la clasif de 4 a 3
    #        if (Table[k,6]==4):
                c[i,j] = 3
                d[i,j] = np.nan
            if Table[k,6]==-10:
                c[i,j] = np.nan
            if Table[k,6]== 4:
                c[i,j] = np.nan
            if (Table[k,6]== 3) and ((Table[k,2] > 20 )and(Table[k,2] < 50)):
                c[i,j] = np.nan
            if Table[k,6]!=4:
                d[i,j] = np.nan
            else:
                d[i,j] = Table[k,7]
            if (Table[k,6]== 3) and ((Table[k,2] > 20 )and(Table[k,2] < 50)):
                e[i,j] = 5
            else:
                e[i,j] = np.nan
            k = k+1
        #c = np.transpose(c)
        #d = np.transpose(d)
        #e = np.transpose(e)
    v1a = masked_array(c, mask=np.isnan(c))
    v1b = masked_array(d, mask=np.isnan(d))
    v1c = masked_array(e, mask=np.isnan(e))
    return v1a,v1b,v1c



def lyap_filter(Table):
    Lyap=Table[:,-1]
    Freq=Table[:,2]
    ghrange=np.array([np.min(Table[:,0]),np.max(Table[:,0])])*1000
    gsdrange=np.array([np.min(Table[:,1]),np.max(Table[:,1])])*1000
    ghN=len(np.unique(Table[:,0]))
    gsdN=len(np.unique(Table[:,1]))
    ratio=np.diff(gsdrange)/np.diff(ghrange)
    Lyap2=np.reshape(Lyap,(ghN,gsdN))
    Freq2=np.reshape(Freq,(ghN,gsdN))
    minFR=2
    maxFR=4.5
    Mask1=0.6 + 0.4* ((Freq2>minFR)*(Freq2<maxFR))
##Colormap functions:
#'''Colormap.set_under(color) Set low out-of-range color.
#Colormap.set_over(color) Set high out-of-range color.
#Colormap.set_bad(color) Set masked-data colo''
    cmap=plt.get_cmap('jet')
    cmap.set_under('k',1)

    FreqN=np.minimum(Freq2,20*np.ones_like(Freq2))/20.  #Maxima freq será 20.
    FreqN[FreqN==0]=-10
    ImgFreq=cmap(FreqN)
    ImgFreq2=ImgFreq*Mask1[:,:,None]

#ImgFreq3=ImgFreq*nonchaosMask[:,:,None]
    ImgLyap=cmap(Lyap2)
    ImgLyap2=ImgLyap*Mask1[:,:,None]
    ImgLyap2[Lyap2==-10,3]=1
    extent = (gsdrange[0],gsdrange[1],ghrange[0],ghrange[1])
    output = [ImgFreq, ImgLyap, ImgLyap2, extent, ratio]
    return output

def Complex_filter(Table):
#    Complex=((Table[:,14]+Table[:,15]+Table[:,16]+Table[:,17])/4.0).tolist()
    Complex=(Table[:,17]).tolist()

    gXrange=np.array([np.min(Table[:,0]),np.max(Table[:,0])])*1000
    gYrange=np.array([np.min(Table[:,1]),np.max(Table[:,1])])*1000
    gXN=len(np.unique(Table[:,0]))
    gYN=len(np.unique(Table[:,1]))

    ratio=np.diff(gYrange)/np.diff(gXrange)

    Complex2=np.reshape(Complex,(gXN,gYN))

    cmap=plt.get_cmap('jet')
    cmap.set_under('k',1.)

    ImgComplex=cmap(Complex2)
    extent = (gYrange[0],gYrange[1],gXrange[0],gXrange[1])
    output = [ImgComplex, extent, ratio]
    return output

def mle_filter(Table):
    if Table.shape[1] == 3:
        lya0=Table[:,2]
        ghrange=np.array([np.min(Table[:,1]),np.max(Table[:,1])])
        gsdrange=np.array([np.min(Table[:,0]),np.max(Table[:,0])])
        ghN=len(np.unique(Table[:,1]))
        gsdN=len(np.unique(Table[:,0]))
        ratio=np.diff(gsdrange)/np.diff(ghrange)
        lya00=np.reshape(lya0,(gsdN,ghN)).transpose()
        extent = (gsdrange[0],gsdrange[1],ghrange[0],ghrange[1])
        output=[lya00,extent]
    else:
        lya0=Table[:,3]
        Freq=Table[:,2]
        gsdrange=np.array([np.min(Table[:,1]),np.max(Table[:,1])])
        ghrange=np.array([np.min(Table[:,0]),np.max(Table[:,0])])
        ghN=len(np.unique(Table[:,1]))
        gsdN=len(np.unique(Table[:,0]))
        ratio=np.diff(gsdrange)/np.diff(ghrange)
        lya00=np.reshape(lya0,(gsdN,ghN))
        Freq00=np.reshape(Freq,(gsdN,ghN))
        extent = (gsdrange[0],gsdrange[1],ghrange[0],ghrange[1])
        output=[Freq00,lya00,ratio,extent]
    return output


def FR_filter(Table):
    Lyap=Table[:,-1]
    Freq=Table[:,2]
    ghrange=np.array([np.min(Table[:,0]),np.max(Table[:,0])])*1000
    gsdrange=np.array([np.min(Table[:,1]),np.max(Table[:,1])])*1000
    ghN=len(np.unique(Table[:,0]))
    gsdN=len(np.unique(Table[:,1]))
    extent=(gsdrange[0],gsdrange[1],ghrange[0],ghrange[1])
    Lyap2=np.reshape(Lyap,(ghN,gsdN))
    Freq2=np.reshape(Freq,(ghN,gsdN))
    output=[Freq2,extent]
    return output

def mask_function(aux0,Mask):
    Freq2=aux0[0]
    Lyap2=aux0[1]
    ratio=aux0[2]
    extent=aux0[3]
    cmap=plt.get_cmap('jet')
    FreqN=np.minimum(Freq2,20*np.ones_like(Freq2))/20.
    ImgFreq=cmap(FreqN)
    ImgFreq2=ImgFreq*Mask[:,:,None]

    ImgLyap=cmap(Lyap2)
    ImgLyap2=ImgLyap*Mask[:,:,None]
    output=[ImgFreq, ImgFreq2,ImgLyap,ImgLyap2, ratio, extent]
    return output


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


def frac_burst(sp_chao,sp_nonchao,N):
    def cyclefunction(sbe):
        cycle=[]
        for i in range(50):
            sp=sbe[i]
            BTime,CT,ev=burst_detect(sp)
    #        print(ev)
            cycle.append(ev)
        return cycle
    
    
    # how to save the data into  two-dimensional (2D) lists in Python
    elem=[]
    elemb=[]
    elements=[]
    elementsb=[]
    for i in range(N):
        elem.append([])
        elemb.append([])
        elements.append([])
        elementsb.append([])
    
    #calculation the burster for all the neurons with different coupling strength g
    for i in range(N):
        for j in range(24):
            p=sp_chao[i][j]
            pp=sp_nonchao[i][j]
            a=cyclefunction(p)  # call the function of finding the number of burst
            aa=cyclefunction(pp)
            elements[i].append(a)
            elementsb[i].append(aa)
    
    for i in range(N):
        for j in range(24):
            per0=[]
            per1=[]
            for k in range(50):
                bbb=elements[i][j][k]
                b4=elementsb[i][j][k]
                per_chaos=float(np.sum(1 for i0 in bbb if i0>1))/len(bbb)
                per_nonchaos=float(np.sum(1 for i0 in b4 if i0>1))/len(b4)
                per0.append(per_chaos)
                per1.append(per_nonchaos)
            elem[i].append(np.mean(per0))
            elemb[i].append(np.mean(per1))
    return elem, elemb