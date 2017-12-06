# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 10:58:52 2016
keshengxuu@gmail.com
@author: ksxuu  
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib 
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D
#The module of generating a colormap index base on discrete intervals
import matplotlib.colors as colors
from PIL import Image


execfile('lib/Defaultsets.py')
plt.rcParams['mathtext.sf'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
#changing the xticks and  yticks fontsize for all sunplots
plt.rc('xtick', labelsize=9) 
plt.rc('ytick', labelsize=9) 
plt.rc('font',size=11)

##################################################
# The classes and functions used for the figures
##################################################
class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def smallworld(Nnode = None, P = None, R = None):
    
    if Nnode is None and P is None and R is None:
        Nnode=50
        p=0.4
        R=0.9

    np.random.seed(10)
    
    angles=np.linspace(0,2*np.pi,num=Nnode,endpoint=False)
    xcoords=np.cos(angles)*R
    ycoords=np.sin(angles)*R
    
    colors=np.random.uniform(0.0,0.3,size=Nnode)
    
    CM = np.zeros((Nnode,Nnode))
    
    for i in range(1,Nnode):
        c_p = np.random.uniform(size=i)
        if c_p[-1]<p:
            CM[i,np.argmax(c_p)]=1
        
        CM[i,i-1]=1
    CM[-1,0]=1
    
    CM = CM+CM.T
    CMl = np.tril(CM)
    
    cc=np.array(np.where(CMl==1)).T
    return xcoords,ycoords,cc, colors

###--------------------------------------------------------------------------#
savepath = 'results/'
Xdata=0
Ydata=4
Ydata2=3
Ggjvals=np.logspace(-6,0,num=23)
Ggjvals=np.insert(Ggjvals,0,10**(-7))
#%%
#%%
def FigX_plots(Burst_chaos,Burst_nonch):
    fig = plt.figure(1, figsize=(8,6))
#    filled_markers = ['o','^']

#%%
#    location=[4,2,4,2]
#  plot SBEs count with different coupling strength
    plt.plot(Ggjvals,Burst_nonch[0],'g-s',markersize=8)
    plt.plot(Ggjvals,Burst_chaos[0],'r--*',  markersize=8)
    plt.legend(('Non-chaotic networks ','Chaotic networks'),loc=2,frameon=False,fontsize=10)
    plt.yticks(np.linspace(0,1,5))
    plt.xscale('log')
    plt.ylabel('Mean fraction of  bursting events',fontsize = 'x-large')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel(r'g', fontsize = 'x-large')
    plt.grid()
    fig.subplots_adjust(left=0.1,bottom=0.1,top=0.985,right=0.985)
    
#    currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
#    parentdir = os.path.dirname(currentdir)
#    sys.path.insert(0, parentdir)
#    os.chdir(parentdir)

    plt.savefig(savepath+'FigX_FraBurst.eps')
    plt.savefig(savepath+'FigX_FraBurst.pdf')
    plt.savefig(savepath+'FigX_FraBurst.png')
    return None


def Fig1_Schematic(times, vt_chaos, vt_nonchaos):
    """
    parameter
    
    times: simulation times
    vt_chaos:  time series(action potentional) of isolate chaotic neuron;
    vt_nonchaos : time series(action potentional) of isolate chaotic neuron;
    """
    #calculation ISIS
    #chaotic
    spikes=np.where(np.diff(1*(vt_chaos[:,0]>-30))>0)[0]*dt
    ISI_chaos=np.diff(spikes)
    #nonchaotic 
    spikes1=np.where(np.diff(1*(vt_nonchaos[:,0]>-30))>0)[0]*dt
    ISI_nonchaos=np.diff(spikes1)
    time1= np.arange(0,100,0.01)
    y1=time1**2
    y2=(time1+5)**2
    
    
    gsd0=0.2846;gh0=0.1829;MLE0=-5.145*10**(-5)
    gsd1=0.203 ;gh1=0.335;MLE1=0.0051

    bb=30
    plt.figure(1, figsize=(12,6))
    plt.clf
    gs = gridspec.GridSpec(100, 100, wspace=0, hspace=0.1)
    ax1=plt.subplot(gs[:19,bb:98])
    ax2=plt.subplot(gs[23:42,bb:98])
    ax3=plt.subplot(gs[57:76,bb:98])
    ax4=plt.subplot(gs[80:99,bb:98])
    ax5=plt.subplot(gs[58:98,0:20])
    ax6=plt.subplot(gs[0:45,0:20])
    
    # examples of chaotic and nonchaotic oscillators
    ax1.plot(times,vt_nonchaos[:,0])
    ax1.set_title(r'Non-chaotic oscillation ($\mathsf{g_{sd}=%0.4f,g_{h}=%0.4f,MLE=%0.4f}$)'%(gsd0,gh0,MLE0),fontsize=13)
    ax1.set_yticks(np.linspace(-90,30,3))
    ax1.set_ylabel(u'V (mv)',fontsize='large')
    ax1.text(-3000,37,'C',fontsize='xx-large')
    ax1.set_xticks([])

    #subplots 2 for ISI of  nonchaos
    ax2.plot(spikes1[1:]/1000,ISI_nonchaos,'.',ms=4)
    ax2.set_xlabel(u'Time (s)',fontsize='large')
    ax2.set_ylabel(u'ISI (ms)',fontsize='large')
    ax2.set_yscale('log')
    ax2.set_ylim([10**2,10**3])
    ax2.set_yticks([])
    ax2.set_yticks([10**2,10**3],['10^2','10^3'])

    #subplots 3 for chaos
    ax3.plot(times,vt_chaos[:,0])
    ax3.set_title(r'Chaotic oscillation ($\mathsf{g_{sd}=%0.4f,g_{h}=%0.4f,MLE=%0.4f}$)'%(gsd1,gh1,MLE1),fontsize=13)
    ax3.set_yticks(np.linspace(-90,30,3))
    ax3.set_ylabel(u'V (mv)',fontsize='large')
    ax3.set_xticks([])

    #subplots 4 for  ISI of chaos
    ax4.plot(spikes[1:]/1000,ISI_chaos,'.',ms=4)
    ax4.set_xlabel(u'Time (s)',fontsize='large')
    ax4.set_ylabel(u'ISI (ms)',fontsize='large')
    ax4.set_yscale('log')
    ax4.set_ylim([10**2,10**3.4])
    ax4.set_yticks([])
    ax4.set_yticks([10**2,10**3],['10^2','10^3'])
    
    #plot the diverage betweebn two trajectionary
    ax5.plot(time1,y1)
    ax5.plot(time1,y2)
    ax5.set_ylim([0,300])
    ax5.set_xlim([5,10])
    ax5.set_axis_off()
    ax5.annotate("", xy=(time1[600],y1[600]),xytext=(time1[610], y2[610]),
        arrowprops=dict(arrowstyle="<-")) 
    ax5.annotate("Neighboring trajectory ", xy=(8,170),xytext=(7, 300),
        arrowprops=dict(arrowstyle="->")) 
    ax5.annotate("", xy=(time1[950],y1[950]),xytext=(time1[960], y2[960]),
        arrowprops=dict(arrowstyle="<-"))
    # The text setting of subplot 5
    ax5.text(time1[600]-0.2,y1[600]-25,'$X_{0}$')
    ax5.text(time1[610]-0.5,y2[610]+25,'$X_{0}+\Delta x_{0} $')
    ax5.text(7,50,' Trajectory',rotation=15)
    ax5.text(time1[950]+0.1,y1[950]+60,' $\Delta x(X_{0},t)$')
    ax5.text(6,-50,'Lyapunov Exponent',fontsize='x-large',color='blue')
    ax5.text(5,300,'B',fontsize='xx-large')
    
    # subplot of small world plot
    # call function of small world
    xcoords,ycoords,cc, colors = smallworld()
    ax6.scatter(xcoords,ycoords,s=80,marker='o',cmap='jet',c=colors,
                linewidths=0,vmin=0,vmax=1,alpha=1,zorder=3)
