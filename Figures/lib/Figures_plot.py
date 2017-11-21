# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 10:58:52 2016
keshengxuu@gmail.com
@author: ksxuu  
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.markers as Mar
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D
#The module of generating a colormap index base on discrete intervals
import matplotlib.colors as colors
execfile('lib/Defaultsets.py')
plt.rcParams['mathtext.sf'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
#changing the xticks and  yticks fontsize for all sunplots
plt.rc('xtick', labelsize=11) 
plt.rc('ytick', labelsize=11) 
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
    for c1,c2 in cc:
        plt.plot((xcoords[c1],xcoords[c2]),(ycoords[c1],ycoords[c2]),'-',
                 color='k')
    ax6.set_xlim((-1,1))
    ax6.set_ylim((-1,1))
    ax6.set_axis_off()
    ax6.text(-0.4,-1.15,'Small World',fontsize='large',color='blue')
    ax6.text(-1.0,1.0,'A',fontsize='xx-large')

    plt.subplots_adjust(bottom=0.08,left=0.02,wspace = 0.4,hspace = 0.25,right=0.97, top=0.95)
    plt.savefig(savepath+'Fig1_Schematic figure.pdf')
    plt.savefig(savepath+'Fig1_Schematic figure.eps')
    return None

def Fig1_regluar(times, vt_chaos, vt_nonchaos):
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
    
    
    gsd0=0.280;gh0=0.181;MLE0=-5.65*10**(-5)
    gsd1=0.203 ;gh1=0.335;MLE1=0.005

    bb=30
    plt.figure(1, figsize=(12,6))
    plt.clf
    gs = gridspec.GridSpec(100, 100, wspace=0, hspace=0.1)
#    ax1=plt.subplot(gs[:19,bb:98])
#    ax2=plt.subplot(gs[23:42,bb:98])
#    ax3=plt.subplot(gs[57:76,bb:98])
#    ax4=plt.subplot(gs[80:99,bb:98])
#    ax5=plt.subplot(gs[58:98,0:20])
#    ax6=plt.subplot(gs[0:45,0:20])

    ax1=plt.subplot(gs[:21,0:71])
    ax2=plt.subplot(gs[22:42,0:71])
    ax3=plt.subplot(gs[57:78,0:71])
    ax4=plt.subplot(gs[79:99,0:71])
    ax5=plt.subplot(gs[54:98,72:98],projection='3d')
    ax6=plt.subplot(gs[0:49,72:98],projection='3d')
    
    # examples of chaotic and nonchaotic oscillators
    ax1.plot(times,vt_nonchaos[:,0])
    ax1.set_title(r'Non-chaotic oscillation ($\mathsf{g_{sd}=%0.3f,g_{h}=%0.3f,MLE=%0.2e}$)'%(gsd0,gh0,MLE0),fontsize='x-large')
    ax1.set_yticks(np.linspace(-90,30,3))
    ax1.set_ylabel(u'V (mV)',fontsize='x-large')
    ax1.text(-800,32,'A',fontsize='xx-large')
    ax1.set_xlim([0,10000])
    ax1.set_xticks([])
    plt.margins(0,0)

    #subplots 2 for ISI of  nonchaos
    ax2.plot(spikes1[1:]/1000,ISI_nonchaos,'.',ms=4)
    ax2.set_xlim([0,10])
    ax2.set_xlabel(u'Time (s)',fontsize='x-large')
    ax2.set_ylabel(u'ISI (ms)',fontsize='x-large')
    ax2.set_yscale('log')
    ax2.set_ylim([10**1.9,10**3.1])
    ax2.set_yticks(np.logspace(2, 3, num=2))
#    ax2.set_yticks([])
#    ax2.set_yticks([100,1000],[r'$10^{2}$',r'$10^{3}$'])

#    ax2.tick_params(axis='y',labelsize=4)
#    ax2.set_yticklabels([r'$10^{2}$',r'$10^{3}$'],  fontsize=8)
    plt.margins(0,0)

    #subplots 3 for chaos
    ax3.plot(times,vt_chaos[:,0])
    ax3.set_title(r'Chaotic oscillation ($\mathsf{g_{sd}=%0.3f,g_{h}=%0.3f,MLE=%0.3f}$)'%(gsd1,gh1,MLE1),fontsize='x-large')
    ax3.set_xlim([0,10000])
    ax3.set_yticks(np.linspace(-90,30,3))
    ax3.set_ylabel(u'V (mV)',fontsize='x-large')
    ax3.text(-800,32,'B',fontsize='xx-large')
    ax3.set_xticks([])
    plt.margins(0,0)

    #subplots 4 for  ISI of chaos
    ax4.plot(spikes[1:]/1000,ISI_chaos,'.',ms=4)
    ax4.set_xlim([0,10])
    ax4.set_xlabel(u'Time (s)',fontsize='x-large')
    ax4.set_ylabel(u'ISI (ms)',fontsize='x-large')
    ax4.set_yscale('log')
    ax4.set_ylim([10**1.9,10**3.1])
    ax4.set_yticks(np.logspace(2, 3, num=2))
    plt.margins(0,0)
#    plt.subplots_adjust(bottom=0.08,left=0.1,wspace = 0.1,hspace = 0.25,right=0.97, top=0.95)
    r'$\mathsf{g_{h}\; (mS/cm^2)}}$'
    # THE 3d PLOTS for regular and chaotic respectively!
    ax5.plot(vt_chaos[6000:,2],vt_chaos[6000:,3],vt_chaos[6000:,4])
    ax5.dist = 8.5
#    ax5.set_title('Chaotic')
    ax6.plot(vt_nonchaos[20000:,2],vt_nonchaos[20000:,3],vt_nonchaos[20000:,4])
    ax6.dist = 8.5

#    ax6.set_title('Regular')

    for ax0 in [ax5,ax6]:
        ax0.set_xlabel(r'$\mathsf{a_{sd}}}$',fontsize='x-large',labelpad=-3)
        ax0.set_ylabel(r'$\mathsf{a_{sr}}}$',fontsize='x-large',labelpad=0)
        ax0.set_zlabel(r'$\mathsf{a_{h}}}$',fontsize='x-large',labelpad=0)
        ax0.tick_params(axis='x',pad=-3)
        ax0.tick_params(axis='y',pad=-1)
        ax0.set_xlim([0.0, 0.65])
        ax0.set_xticks([0.0, 0.2, 0.4, 0.6])
        ax0.set_xticklabels(['0.0', '0.2', '0.4', '0.6'])
        
        ax0.set_yticks([0.2, 0.3, 0.4])
        ax0.set_yticklabels(['0.2', '0.3', '0.4'])


    ax5.set_ylim([0.19,0.4])
    ax6.set_ylim([0.15,0.4])
    ax5.set_zlim([0.025,0.12])
    ax5.set_zticks([0.03, 0.06, 0.09, 0.12])
    ax5.set_zticklabels([ '0.03', '0.06', '0.09','0.12'])


    ax6.set_zlim([0.025,0.16])
    ax6.set_zticks([0.04, 0.08, 0.12, 0.16])
    ax6.set_zticklabels([ '0.04', '0.08', '0.12','0.16'])


    plt.subplots_adjust(bottom=0.08,left=0.06,wspace = 0.1,hspace = 0.1,right=0.97, top=0.96)
    plt.savefig(savepath+'Fig1_regu.pdf')
    plt.savefig(savepath+'Fig1_regu.png')
    plt.savefig(savepath+'Fig1_regu.eps')
    return None


def  Fig2_FS(mle, FR, FR_mask,tupTable2):
    """
    parameter
    
    MLE : maximum lyapuonv expoents
    FR : firing rate
    FR_MASK : tuple
    
    tupTable2: tuple
         should contain
          chaor1,nonchaor1, chaor2, nonchaor2 

    """
#    Table, Table_FP, Table_FR = tupTable 
    chaor1,nonchaor1, chaor2, nonchaor2 = tupTable2 
    v1a,v1b,v1c = FR_mask
    
    chaos=[chaor1,chaor2]
    nonchaos=[nonchaor1,nonchaor2]

    cmap=plt.get_cmap('jet')
    cmap.set_under('k',1)
    norm=colors.Normalize(vmin=0,vmax=0.006)
    fig=plt.figure(1,figsize=(10,7))
    plt.clf()
    gs = gridspec.GridSpec(60, 90, wspace=0, hspace=0.1)
    ax1=plt.subplot(gs[0:20,0:22])
    ax2=plt.subplot(gs[0:20,32:54])
    ax3=plt.subplot(gs[0:20,65:87])
    ax4=plt.subplot(gs[27:42,0:24])
    ax5=plt.subplot(gs[27:42,32:56])
    ax6=plt.subplot(gs[27:42,65:89])
    ax7=plt.subplot(gs[45:60,0:24])
    ax8=plt.subplot(gs[45:60,32:56])
    ax9=plt.subplot(gs[45:60,65:89])
    
    ax_iso = [ax1,ax2,ax3]
    ax = [ax4,ax5,ax6,ax7,ax8,ax9]
    #subplots of MLE
    im=ax1.imshow(mle[0],origin='lower',extent=mle[1],norm=norm,aspect='auto',interpolation='none',cmap=cmap)
    ax1.autoscale(False)
    ax1.plot((gsdmin1,gsdmin1,gsdmax1,gsdmax1,gsdmin1),(ghmin1,ghmax1,ghmax1,ghmin1,ghmin1),'r-',lw=2)
    ax1.plot((gsdmin2,gsdmin2,gsdmax2,gsdmax2,gsdmin2),(ghmin2,ghmax2,ghmax2,ghmin2,ghmin2),'g-',lw=2)
    ax1.plot((gsdmin1b,gsdmin1b,gsdmax1b,gsdmax1b,gsdmin1b),(ghmin1b,ghmax1b,ghmax1b,ghmin1b,ghmin1b),'r-',lw=2)
    ax1.plot((gsdmin2b,gsdmin2b,gsdmax2b,gsdmax2b,gsdmin2b),(ghmin2b,ghmax2b,ghmax2b,ghmin2b,ghmin2b),'g-',lw=2)
    ax1.text(0.11,0.58,'A',fontsize='xx-large')
    #the colorbar for MLE
    cax1=fig.add_axes([0.32, 0.68, 0.01, 0.26])
    cbar1=fig.colorbar(im, extend='both',cax=cax1)
    cbar1.set_label(r'MLE ',fontsize='large',labelpad=-20, y=1.1, rotation=0)
    cbar1.set_ticks(np.linspace(0,0.006,4))
    cbar1.ax.tick_params(labelsize='x-small')

    #subplots of firing rate
    im2=ax2.imshow(FR[0],origin='lower',extent=FR[1],
               aspect='auto',vmin=0,vmax=20,interpolation='none',cmap=cmap)
    ax2.autoscale(False) # to reset the axes scale.
    ax2.plot((gsdmin1,gsdmin1,gsdmax1,gsdmax1,gsdmin1),(ghmin1,ghmax1,ghmax1,ghmin1,ghmin1),'r-',lw=2)
    ax2.plot((gsdmin2,gsdmin2,gsdmax2,gsdmax2,gsdmin2),(ghmin2,ghmax2,ghmax2,ghmin2,ghmin2),'g-',lw=2)
    ax2.plot((gsdmin1b,gsdmin1b,gsdmax1b,gsdmax1b,gsdmin1b),(ghmin1b,ghmax1b,ghmax1b,ghmin1b,ghmin1b),'r-',lw=2)
    ax2.plot((gsdmin2b,gsdmin2b,gsdmax2b,gsdmax2b,gsdmin2b),(ghmin2b,ghmax2b,ghmax2b,ghmin2b,ghmin2b),'g-',lw=2)
    # the setting of firing rate colorbar
    cax2=fig.add_axes([0.64, 0.68, 0.01, 0.26])
    cbar2=fig.colorbar(im2, extend='max',cax=cax2)
    cbar2.set_label('FR ',fontsize='large',labelpad=-18, y=1.1, rotation=0)
    cbar2.set_ticks(np.linspace(0,20,5))
    cbar2.ax.tick_params(labelsize='small')
    
    # the subplots of Firing pattern
    cmaps = mpl.colors.ListedColormap(['black','blue','palegreen','yellow','darkorange','brown','darkred'])
    im3a= ax3.imshow(v1b, origin = 'lower', vmin=1, vmax=10, extent=mle[1],
                     interpolation='nearest', cmap=cm.Oranges, norm = MidpointNormalize(midpoint=2), aspect='auto')
    im3b = ax3.imshow(v1a, origin = 'lower', vmin=0, vmax=6, extent=mle[1],
                      interpolation='nearest',cmap = cmaps,  aspect='auto')#, alpha = 0.5)
    ax3.autoscale(False) # to reset the axes scale.
    ax3.plot((gsdmin1,gsdmin1,gsdmax1,gsdmax1,gsdmin1),(ghmin1,ghmax1,ghmax1,ghmin1,ghmin1),'r-',lw=2)
    ax3.plot((gsdmin2,gsdmin2,gsdmax2,gsdmax2,gsdmin2),(ghmin2,ghmax2,ghmax2,ghmin2,ghmin2),'g-',lw=2)
    ax3.plot((gsdmin1b,gsdmin1b,gsdmax1b,gsdmax1b,gsdmin1b),(ghmin1b,ghmax1b,ghmax1b,ghmin1b,ghmin1b),'r-',lw=2)
    ax3.plot((gsdmin2b,gsdmin2b,gsdmax2b,gsdmax2b,gsdmin2b),(ghmin2b,ghmax2b,ghmax2b,ghmin2b,ghmin2b),'g-',lw=2)
    # the colorbar of the firing pattern
    cax3 = fig.add_axes([0.965, 0.68, 0.012, 0.26])
    cbar3 = fig.colorbar(im3b, ticks=[ 0.4, 1.3, 2.1, 3, 3.8, 4.7, 5.5],cax=cax3)
    cbar3.ax.set_yticklabels(['0', '1', '2', '3', '4', '5', '6'])
    cbar3.set_label(r'FP',fontsize='large',labelpad=-18, y=1.1, rotation=0)
    cbar3.ax.tick_params(labelsize='medium')
    
    # the plots for synchronization transition,inlcude metastablilty and MLE 
    #subplots of the order parameters for the range1 and range2 repectively.
    for i,j in zip(range(0,4,3),range(0,2)):
        ax[i].plot(nonchaos[j][:,0][1],np.mean(nonchaos[j][:,1:11],axis=1)[1],'g-s', markersize=4.5)
        ax[i].plot(nonchaos[j][:,0][5:],np.mean(nonchaos[j][:,1:11],axis=1)[5:],'g-s', markersize=4.5,
          label='Non-chaotic networks ')
        
        ax[i].plot(chaos[j][:,0][1],np.mean(chaos[j][:,1:11],axis=1)[1],'r--*',  markersize=5)
        ax[i].plot(chaos[j][:,0][5:],np.mean(chaos[j][:,1:11],axis=1)[5:],'r--*',  markersize=5,
          label='Chaotic networks')
        ax[i].set_xscale('log')
        ax[i].set_ylabel(r'R', fontsize = 'large')
        ax[i].set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
        ax[i].set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'])
        y1,y2=ax[i].get_ylim()
        ax[i].bar(3e-6,0.5,3e-6,-0.01,color='white',lw=5,clip_on=False,zorder=5)
        ax[i].set_ylim((y1,y2))

        ax[i].set_xlim((5e-7,1.5))
        ax[i].set_xticks((1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1))
        ax[i].set_xticklabels(("$0$","$10^{-5}$","$10^{-4}$","$10^{-3}$","$10^{-2}$","$10^{-1}$","$10^0$"))
        
        t1=ax[i].xaxis.get_major_ticks()[0]
        t1.set_pad(5.5)
    
        for t in ax[i].xaxis.get_minor_ticks()[:9]:
            t.set_visible(False)
        ax[i].tick_params(labelsize=10)
        

    ax4.legend(loc=2,frameon=False,fontsize=8)
        
    ax[0].text(10**(-8.5),1.1,'B',fontsize='xx-large')
#    ax[3].text(10**(-10),1.1,'C',fontsize='xx-large')


    #subplots of mestablilty
    for i,j in zip(range(1,5,3),range(0,2)):
        ax[i].plot(nonchaos[j][:,0][1],np.mean(nonchaos[j][:,21:31],axis=1)[1],'g-s', markersize=4.5)
        ax[i].plot(nonchaos[j][:,0][5:],np.mean(nonchaos[j][:,21:31],axis=1)[5:],'g-s', markersize=4.5)

        ax[i].plot(chaos[j][:,0][1],np.mean(chaos[j][:,21:31],axis=1)[1],'r--*',  markersize=5)
        ax[i].plot(chaos[j][:,0][5:],np.mean(chaos[j][:,21:31],axis=1)[5:],'r--*',  markersize=5)

        ax[i].set_ylabel(u'Metastability', fontsize='large',labelpad= -1)
        ax[i].set_ylim([-0.001,0.045])
        ax[i].set_yticks(np.linspace(0.000,0.045,4))
        ax[i].set_xscale('log')
        y1,y2=ax[i].get_ylim()
        ax[i].bar(3e-6,0.1,3e-6,-0.1,color='white',lw=5,clip_on=False,zorder=5)
        ax[i].set_ylim((y1,y2))

        ax[i].set_xlim((5e-7,1.5))
        ax[i].set_xticks((1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1))
        ax[i].set_xticklabels(("$0$","$10^{-5}$","$10^{-4}$","$10^{-3}$","$10^{-2}$","$10^{-1}$","$10^0$"))
        
        t1=ax[i].xaxis.get_major_ticks()[0]
        t1.set_pad(5.5)
    
        for t in ax[i].xaxis.get_minor_ticks()[:9]:
            t.set_visible(False)
        ax[i].tick_params(labelsize=10)



    
    
    
    #subplots for the the mle
    Ran=['Region 1','Region 2']
    for i,j,Ra in zip(range(2,6,3),range(0,2),Ran):
        ax[i].plot(nonchaos[j][:,0][1],np.mean(nonchaos[j][:,11:21],axis=1)[1],'g-s', markersize=4.5, color='g')
        ax[i].plot(nonchaos[j][:,0][5:],np.mean(nonchaos[j][:,11:21],axis=1)[5:],'g-s', markersize=4.5, color='g')
        
        ax[i].plot(chaos[j][:,0][1],np.mean(chaos[j][:,11:21],axis=1)[1],'r--*',  markersize=5)
        ax[i].plot(chaos[j][:,0][5:],np.mean(chaos[j][:,11:21],axis=1)[5:],'r--*',  markersize=5)
        ax[i].set_ylabel(u'MLE', fontsize='large')
        ax[i].set_ylim([-0.0005,0.0122])
        ax[i].set_yticks(np.linspace(0.000,0.012,4))
        ax[i].text(10**(0.35),0.007,' %s'% Ra,rotation=90,fontsize='large')
        ax[i].set_xscale('log')
        y1,y2=ax[i].get_ylim()
        ax[i].bar(3e-6,0.1,3e-6,-0.1,color='white',lw=5,clip_on=False,zorder=5)
        ax[i].set_ylim((y1,y2))

        ax[i].set_xlim((5e-7,1.5))
        ax[i].set_xticks((1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1))
        ax[i].set_xticklabels(("$0$","$10^{-5}$","$10^{-4}$","$10^{-3}$","$10^{-2}$","$10^{-1}$","$10^0$"))
        
        t1=ax[i].xaxis.get_major_ticks()[0]
        t1.set_pad(5.5)
    
        for t in ax[i].xaxis.get_minor_ticks()[:9]:
            t.set_visible(False)
        ax[i].tick_params(labelsize=10)


        
    
    # the setting for text and xticks and y ticks FOR subplots 1-3
    ax1.set_ylabel(r'$\mathsf{g_{h}\; (mS/cm^2)}}$', fontsize = 'large')
    for ax0 in ax_iso:
            ax0.text(0.208,0.36,r'1',color='r',fontsize='large',fontweight='bold')
            ax0.text(0.208,0.11,r'1',color='g',fontsize='large',fontweight='bold')
            ax0.text(0.29,0.36,r'2',color='r',fontsize='large',fontweight='bold')
            ax0.text(0.29,0.19,r'2',color='g',fontsize='large',fontweight='bold')
            ax0.set_xlabel(r'$\mathsf{g_{sd} \; (mS/cm^2)}}$', fontsize = 'large')
            ax0.set_xticks(np.linspace(0.17,0.33,5))
            ax0.set_yticks(np.linspace(0,0.6,4))

    for i in range(3,6):
        ax[i].set_xlabel(r'g', fontsize = 'large')


    plt.subplots_adjust(left = 0.1,bottom=0.1, right=0.98, top=0.97, wspace=0.35, hspace=0.25)
    plt.show
    plt.savefig(savepath+'Fig2_ghsdFS.eps')
    plt.savefig(savepath+'Fig2_ghsdFS.png')
    plt.savefig(savepath+'Fig2_ghsdFS.pdf')
    
    return None

def Fig3Histgram_FS(tupTable, Burst_chaos, Burst_nonch):
    """
    parameter
    
    tupTable : tuple
         should contain desired table of fixed size
          chaosTableFR, nonchaosTableFR, chaosTableFR2, nonchaosTableFR2 

    """
    chaosTableFR, nonchaosTableFR, chaosTableFR2, nonchaosTableFR2  = tupTable
    #Firing rate histogram of range 1 and range2
    fig = plt.figure(3,figsize=(11,3.5))
    plt.clf()
    gs = gridspec.GridSpec(1, 100, wspace=0, hspace=0.1)


    ax1= plt.subplot(gs[0,0:30])
    plt.plot(Ggjvals[1],Burst_nonch[0][1],'g-s',markersize=6)
    plt.plot(Ggjvals[5:],Burst_nonch[0][5:],'g-s',markersize=6,label= 'Non-chaotic networks ')
    
    plt.plot(Ggjvals[1],Burst_chaos[0][1],'r--*',  markersize=6)
    plt.plot(Ggjvals[5:],Burst_chaos[0][5:],'r--*',  markersize=6,label= 'Chaotic networks')
    plt.legend(loc=2,frameon=False,fontsize='small')
    plt.yticks(np.linspace(0,1,5))
    plt.xscale('log')
#    plt.ylabel('Mean fraction of  bursting events',fontsize = 'x-large')
    plt.ylabel('MB',fontsize = 'large')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel(r'g', fontsize = 'large',labelpad = -3)
    y1,y2=ax1.get_ylim()
    ax1.bar(3e-6,0.4,3e-6,-0.1,color='white',lw=0,clip_on=False,zorder=3)
    ax1.set_ylim((y1,y2))

    ax1.set_xlim((5e-7,1.99))
    ax1.set_xticks((1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1))
    ax1.set_xticklabels(("$0$","$10^{-5}$","$10^{-4}$","$10^{-3}$","$10^{-2}$","$10^{-1}$","$10^0$",))

    t1=ax1.xaxis.get_major_ticks()[0]
    t1.set_pad(5.5)
    
    for t in ax1.xaxis.get_minor_ticks()[:8]:
        t.set_visible(False)
        
    plt.text(10**(-7.8), 1 ,'A',fontsize='x-large')
#    plt.grid()
    
    
    labels=['Non-chaotic neurons','Chaotic neurons']
    ax2 = plt.subplot(gs[0,40:68])
    ax2.hist((nonchaosTableFR[:,2],chaosTableFR[:,2]),bins=25,color=('g','r'),label=labels)
    ax2.set_ylabel(u'Event count',fontsize='large')
    ax2.legend(fontsize='small',loc=2,frameon=False)
#    ax2.legend(loc='upper center',bbox_to_anchor=(1.08, 1.18),ncol=10,handleheight=1,handlelength=5,fontsize='small',frameon=False)
    ax2.set_yticks(np.linspace(0,36,4))
    plt.text(1.65, 36 ,'B',fontsize='x-large')
    ax2.set_xlabel(u'Firing rate (spikes/s)',fontsize='large')
#    ax1.text(1.8,20,'B',fontsize='xx-large')
    ax2.set_title('Region 1',fontsize='small')
    
    ax3=plt.subplot(gs[0,71:99])
    ax3.hist((chaosTableFR2[:,2],nonchaosTableFR2[:,2]),bins=25,color=('r','g'),label=labels)
    ax3.set_xlabel(u'Firing rate (spikes/s)',fontsize='large')
    ax3.set_yticks(np.linspace(0,36,4))
    ax3.set_title('Region 2',fontsize='small')
    
    plt.subplots_adjust(left = 0.08,bottom=0.2, right=0.98, top=0.945, wspace=0.20, hspace=0.25)
    plt.show
    plt.savefig(savepath+'Fig3_HistFS.eps')
    plt.savefig(savepath+'Fig3_HistFS.png')
    plt.savefig(savepath+'Fig3_HistFS.pdf')
    return None

def Fig4_Dis(lyap, lyap2, auxchao, auxnochaos, auxnoih, tupTable):
    chaosTable, nonchaosTable2, nonIhTable2 = tupTable
    cmap=plt.get_cmap('jet')
#    cmap.set_under('k',1)
    norm=colors.Normalize(vmin=0,vmax=0.006)
    
    fig=plt.figure(1,figsize=(8,7))
    plt.clf()
    gs = gridspec.GridSpec(60, 80,wspace=0.2,hspace=0.35)
    
    #subplot1
    #sizes = ['xx-small', 'x-small', 'small', 'medium', 'large','x-large', 'xx-large']
    ax1=plt.subplot(gs[0:15,1:20])
    im1=ax1.imshow(lyap[1],origin='lower',extent=lyap[3],aspect='auto',norm=norm,interpolation='none',cmap=cmap)
    ax1.set_xlabel(r'$\mathsf{ g_{sr} \; (mS/cm^2)}$',fontsize='large')
    ax1.set_ylabel(r'$\mathsf{ g_{sd} \; (mS/cm^2)}$',fontsize='large')
    plt.text(0.23,0.415,r'MLE  (HB+$\mathsf{ \; I_{h}}$)',fontsize='medium')
    plt.text(0.15,0.4,r'A',fontsize='x-large')
    
    #subplot2
    ax2=plt.subplot(gs[0:15,35:54])
    im2=plt.imshow(auxchao[1],origin='lower',extent=auxchao[-1],aspect='auto',vmin=0,vmax=20,interpolation='none',cmap=cmap)
    plt.xlabel(r'$\mathsf{ g_{sr} \; (mS/cm^2)}$',fontsize='large')
    plt.ylabel(r'$\mathsf{ g_{sd} \; (mS/cm^2)}$',fontsize='large')
    plt.text(0.215,0.415,r'FR - chaotic (HB+$\mathsf{ \; I_{h}}$)',fontsize='medium')
    
    #subplot3
    ax3=plt.subplot(gs[0:15,58:77])
    im3=plt.imshow(auxnochaos[1],origin='lower',extent=auxnochaos[-1],aspect='auto',vmin=0,vmax=20 ,interpolation='none',cmap=cmap)
    plt.xlabel(r'$\mathsf{ g_{sr} \; (mS/cm^2)}$',fontsize='large')
    plt.text(0.20,0.415,r'FR - non-chaotic (HB+$\mathsf{ \; I_{h}}$)',fontsize='medium')
    
    #subplot4
    ax4=plt.subplot(gs[24:39,1:20])
    im4=plt.imshow(lyap2[1],origin='lower',extent=lyap2[3],aspect='auto',norm=norm,interpolation='none',cmap=cmap)
    plt.ylabel(r'$\mathsf{ g_{sd} \; (mS/cm^2)}$',fontsize='large')
    plt.text(0.235,0.41,r'MLE (HB)',fontsize='small')
    plt.text(0.15,0.4,r'B',fontsize='x-large')
    
    ##subplot5
    ax5=plt.subplot(gs[45:60,1:20])
    im5=plt.imshow(auxnoih[1],origin='lower',extent=auxnoih[-1],aspect='auto',vmin=0,vmax=20,interpolation='none',cmap=cmap)
    plt.xlabel(r'$\mathsf{ g_{sr} \; (mS/cm^2)}$',fontsize='large')
    plt.ylabel(r'$\mathsf{ g_{sd} \; (mS/cm^2)}$',fontsize='large')
    plt.text(0.22,0.42,r'FR (without $\mathsf{ \; I_{h}}$)',fontsize='medium')
    
    #subplot6
    ax6=plt.subplot(gs[29:58,35:80])
    labels=['Non-chaotic','Chaotic','NoIh']
    l1,l2,l3=plt.hist((nonchaosTable2[:,2],chaosTable[:,2],nonIhTable2[:,2]),bins=15,color=('r','g','b'),label=labels)
    plt.yticks(np.linspace(0,50,6))
    plt.ylabel(u'Event count ',fontsize='large')
    plt.legend(loc='upper center',bbox_to_anchor=(0.5, 1.1),ncol=3,fontsize='medium',frameon=False)
    plt.xlabel(r'Firing rate (spikes/s)',fontsize='large')
    plt.xticks(np.linspace(3.0,4.5,4))
    plt.text(2.8,53,r'C',fontsize='x-large')
    ax6.tick_params(labelsize=9)
    
    #fig.legend([l1, l2, l3], ['IC','INC','NoIh'], bbox_to_anchor=[0.5, 0.5], loc='center', ncol=2)
    
    #set ticks and labels
    axies=[ax1,ax2,ax3,ax4,ax5]
    
    for ax in axies:
        ax.set_xlim([0.2,0.35])
        ax.set_xticks(np.linspace(0.2,0.35,4))
        ax.set_yticks(np.linspace(0.10,0.40,4))
        ax.tick_params(labelsize=9)
    
    
    cax1 = fig.add_axes([0.31, 0.745, 0.009, 0.2])
    cbar1=fig.colorbar(im1, extend='both',cax=cax1)
    cbar1.set_label('MLE ',fontsize='large',labelpad=-26,y=1.15,rotation=0)
    cbar1.set_ticks(np.linspace(0,0.006,4))
    cbar1.ax.tick_params(labelsize=9)
    
    
    cax4 = fig.add_axes([0.31, 0.385, 0.009, 0.2])
    cbar4=fig.colorbar(im4, extend='both',cax=cax4)
    cbar4.set_label('MLE ',fontsize='large',labelpad=-26,y=1.15,rotation=0)
    cbar4.set_ticks(np.linspace(0,0.006,4))
    cbar4.ax.tick_params(labelsize=9)
    
    cax2 = fig.add_axes([0.945, 0.745, 0.009, 0.2])
    cbar2=fig.colorbar(im2, extend='max',cax=cax2)
    cbar2.set_label('FR ',fontsize='large',labelpad=-17,y=1.11,rotation=0)
    cbar2.set_ticks(np.linspace(0,20,5))
    cbar2.ax.tick_params(labelsize=9)
    
    cax5 = fig.add_axes([0.31,0.08, 0.01, 0.2])
    cbar5=fig.colorbar(im5, extend='max',cax=cax5)
    cbar5.set_label('FR ',fontsize='large',labelpad=-17,y=1.11,rotation=0)
    cbar5.set_ticks(np.linspace(0,20,5))
    cbar5.ax.tick_params(labelsize=9)
    
    #fig.tight_layout()
    plt.subplots_adjust(left = 0.08,bottom=0.08, right=0.97, top=0.965, wspace=0.0, hspace=0)
    plt.savefig(savepath+'Fig4_gsdgsrDis.eps')
    plt.savefig(savepath+'Fig4_gsdgsrDis.pdf')
    plt.savefig(savepath+'Fig4_gsdgsrDis.png')
    plt.show
    #
    #directory='Data_sets'
    #if not os.path.exists(directory):
    #    os.makedirs(directory)
        
    #np.savetxt(directory+'/FR%sto%schaos.txt'%(10*minFR,10*maxFR),chaosvalues,fmt='%.6f', delimiter=' ',header='gsd\t gsr')
    #np.savetxt(directory+'/FR%sto%snonchaos.txt'%(10*minFR,10*maxFR),nonchaosvalues,fmt='%.6f', delimiter=' ',header='gsd\t gsr')
    #np.savetxt(directory+'/FR%sto%snoih.txt'%(10*minFR,10*maxFR),nonIhvalues,fmt='%.6f', delimiter=' ',header='gsd\t gsr')
    return None

def Fig5_plot(chaos, nonchaos, noIh):
    fig=plt.figure(1,figsize=(11, 5))
    plt.clf
    ax=[plt.subplot(2,3,i) for i in range(1,7)]
#    ax = GridSpec(2,3)
    
    #subplots for the order parameter 
    for i,j in zip(range(0,4,3),range(0,2)):
        ax[i].plot(nonchaos[j][:,0][1],np.mean(nonchaos[j][:,1:6],axis=1)[1],'g-s',label=u'Non-chaotic networks')
        ax[i].plot(nonchaos[j][:,0][5:],np.mean(nonchaos[j][:,1:6],axis=1)[5:],'g-s')
        
        ax[i].plot(chaos[j][:,0][1],np.mean(chaos[j][:,1:6],axis=1)[1],'r--*',label=u'Chaotic networks')
        ax[i].plot(chaos[j][:,0][5:],np.mean(chaos[j][:,1:6],axis=1)[5:],'r--*')
        
        ax[i].plot(noIh[j][:,0][1],np.mean(noIh[j][:,1:6],axis=1)[1],'b-.^',label=u'NoIh networks')
        ax[i].plot(noIh[j][:,0][5:],np.mean(noIh[j][:,1:6],axis=1)[5:],'b-.^')
#        ax[i].tick_params(labelsize=10)
#        ax[i].set_ylim([0.000,1.02])
        ax[i].set_ylabel(r'R', fontsize='large')
        ax[i].set_xscale('log')

        
        y1,y2=ax[i].get_ylim()
        ax[i].bar(3e-6,0.5,3e-6,-0.01,color='white',lw=5,clip_on=False,zorder=5)
        ax[i].set_ylim((y1,y2))
        
        ax[i].set_xlim((5e-7,1.5))
        ax[i].set_xticks((1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1))
        ax[i].set_xticklabels(("$0$","$10^{-5}$","$10^{-4}$","$10^{-3}$","$10^{-2}$","$10^{-1}$","$10^0$"))
        
        t1=ax[i].xaxis.get_major_ticks()[0]
        t1.set_pad(5.5)
        
        for t in ax[i].xaxis.get_minor_ticks()[:8]:
            t.set_visible(False)
        ax[i].tick_params(labelsize=10)




    ax[0].legend(fontsize='small',loc=2,frameon=False)
    ax[0].text(10**(-7.7),0.975,'A',fontsize='xx-large')
    ax[3].text(10**(-7.7),0.975,'B',fontsize='xx-large')
           
    #subplots for the KYdim
    for i,j in zip(range(1,5,3),range(0,3)):
        ax[i].plot(nonchaos[j][:,0][1],np.mean(nonchaos[j][:,6:11],axis=1)[1],'g-s',)
        ax[i].plot(nonchaos[j][:,0][5:],np.mean(nonchaos[j][:,6:11],axis=1)[5:],'g-s')
        
        ax[i].plot(chaos[j][:,0][1],np.mean(chaos[j][:,6:11],axis=1)[1],'r--*')
        ax[i].plot(chaos[j][:,0][5:],np.mean(chaos[j][:,6:11],axis=1)[5:],'r--*')

        ax[i].plot(noIh[j][:,0][1],np.mean(noIh[j][:,6:11],axis=1)[1],'b-.^')
        ax[i].plot(noIh[j][:,0][5:],np.mean(noIh[j][:,6:11],axis=1)[5:],'b-.^')

        ax[i].set_ylabel(u'Metastability', fontsize='large',labelpad=5)
        ax[i].set_xscale('log')
        ax[i].set_yticks(np.linspace(0.000,0.030,4))
        
        y1,y2=ax[i].get_ylim()
        ax[i].bar(3e-6,0.11,3e-6,-0.1,color='white',lw=5,clip_on=False,zorder=5)
        ax[i].set_ylim((y1,y2))

        ax[i].set_xlim((5e-7,1.5))
        ax[i].set_xticks((1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1))
        ax[i].set_xticklabels(("$0$","$10^{-5}$","$10^{-4}$","$10^{-3}$","$10^{-2}$","$10^{-1}$","$10^0$"))
        
        t1=ax[i].xaxis.get_major_ticks()[0]
        t1.set_pad(5.5)
    
        for t in ax[i].xaxis.get_minor_ticks()[:9]:
            t.set_visible(False)
        ax[i].tick_params(labelsize=10)
#        ax[i].tick_params(labelsize=10)
    
    
    
    #subplots for the the mle
    Fre=[[3.0,4.5],[7.0,9.5]]
    for i,j,fr in zip(range(2,6,3),range(0,2),Fre):
        ax[i].plot(nonchaos[j][:,0][1],np.mean(nonchaos[j][:,21:26],axis=1)[1],'g-s')
        ax[i].plot(nonchaos[j][:,0][5:],np.mean(nonchaos[j][:,21:26],axis=1)[5:],'g-s')

        ax[i].plot(chaos[j][:,0][1],np.mean(chaos[j][:,21:26],axis=1)[1],'r--*')
        ax[i].plot(chaos[j][:,0][5:],np.mean(chaos[j][:,21:26],axis=1)[5:],'r--*')

        ax[i].plot(noIh[j][:,0][1],np.mean(noIh[j][:,21:26],axis=1)[1],'b-.^')
        ax[i].plot(noIh[j][:,0][5:],np.mean(noIh[j][:,21:26],axis=1)[5:],'b-.^')

        ax[i].set_ylabel(u'MLE', fontsize='large')
        ax[i].set_ylim([-0.0005,0.0108])
        ax[i].set_yticks(np.linspace(0.000,0.012,4))
        ax[i].set_xscale('log')
        ax[i].text(10**(0.4),0.01,'FR = %g to %g spikes/s'%(fr[0],fr[1]),rotation=90,fontsize='medium')
        y1,y2=ax[i].get_ylim()
        ax[i].bar(3e-6,0.1,3e-6,-0.1,color='white',lw=5,clip_on=False,zorder=5)
        ax[i].set_ylim((y1,y2))

        ax[i].set_xlim((5e-7,1.5))
        ax[i].set_xticks((1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1))
        ax[i].set_xticklabels(("$0$","$10^{-5}$","$10^{-4}$","$10^{-3}$","$10^{-2}$","$10^{-1}$","$10^0$"))
        
        t1=ax[i].xaxis.get_major_ticks()[0]
        t1.set_pad(5.5)
    
        for t in ax[i].xaxis.get_minor_ticks()[:9]:
            t.set_visible(False)
        ax[i].tick_params(labelsize=10)




    #set the xlabel    
    for i in range(3,6):
        ax[i].set_xlabel(r'g', fontsize = 'x-large')
        
                       
#    fig.legend([line1, line2,line3], ['Non-chaotic networks','Chaotic networks', 'NoIh networks'], bbox_to_anchor=[0.525, 1.004], 
#              loc='upper center', ncol=4) 
     
    plt.subplots_adjust(bottom=0.12,left=0.06,wspace = 0.3,hspace = 0.15,right=0.975, top=0.98)
    plt.draw()
    plt.savefig(savepath+'Fig5_gsdgsrDisSyn.eps')
    plt.savefig(savepath+'Fig5_gsdgsrDisSyn.png')
    plt.savefig(savepath+'Fig5_gsdgsrDisSyn.pdf')
    return None
