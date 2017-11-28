# -*- coding: utf-8 -*-
"""
Firing rate statistics
FIXED DISTRIBUTION OF gsd, gh

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib

plt.rcParams['mathtext.sf'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')

#load the datatxt
Table=np.loadtxt("data/mleFR_finalTableghgsd.txt",skiprows=1,delimiter=',')
#load the order parameter data
chaor1= np.loadtxt('data/chaoticdata.txt',delimiter=',')
nonchaor1= np.loadtxt('data/nonchaoticdata.txt',delimiter=',')

chaos=np.loadtxt('data/Chaosghsd_node50_table.txt',delimiter=',')
nonchaos=np.loadtxt('data/Nonchghsd_node50_table.txt',delimiter=',')

chao_da1=chaor1[:,1:]
nonchao_da1=nonchaor1[:,1:]


# the mean value over the different seed
chaomean=np.mean(chao_da1,axis=1)
nonchaomean=np.mean(nonchao_da1,axis=1)

minFR=2
maxFR=3.5
upmle=0.001
lowmle=0.0001
sim=50

Lyap=Table[:,-1]
Freq=Table[:,2]
ghrange=np.array([np.min(Table[:,0]),np.max(Table[:,0])])
gsdrange=np.array([np.min(Table[:,1]),np.max(Table[:,1])])
ghN=len(np.unique(Table[:,0]))
gsdN=len(np.unique(Table[:,1]))
ratio=np.diff(gsdrange)/np.diff(ghrange)
Lyap2=np.reshape(Lyap,(ghN,gsdN))
Freq2=np.reshape(Freq,(ghN,gsdN))


extent = (gsdrange[0],gsdrange[1],ghrange[0],ghrange[1])

chaosTable=Table[(Table[:,2]>minFR)*(Table[:,2]<maxFR)*(Table[:,-1]>upmle)]
nonchaosTable=Table[(Table[:,2]>minFR)*(Table[:,2]<maxFR)*(Table[:,-1]<lowmle)]
nonchaosTable2=[]


for freq in chaosTable[:,2]:
    if freq in nonchaosTable[:,2]:
        indices=np.where(nonchaosTable[:,2]==freq)[0]
        nonchaosTable2.append(nonchaosTable[np.random.choice(np.where(nonchaosTable[:,2]==freq)[0],1)[0]])
    else:
        nonchaosTable2.append(nonchaosTable[np.abs(nonchaosTable[:,2]-freq).argmin()])


nonchaosTable2=np.array(nonchaosTable2)

chaosindex=np.random.choice(len(chaosTable),size=sim)
chaosvalues=chaosTable[chaosindex,0:2]

nonchaosindex=np.random.choice(len(nonchaosTable2),size=sim)
nonchaosvalues=nonchaosTable2[nonchaosindex,0:2]


chaosMask=0.5 * np.ones_like(Freq2) + 0.5* ((Freq2>minFR)*(Freq2<maxFR)*(Lyap2>upmle))
nonchaosMask=0.5 * np.ones_like(Freq2) + 0.5* ((Freq2>minFR)*(Freq2<maxFR)*(Lyap2<lowmle))


cmap=plt.get_cmap('jet')
FreqN=np.minimum(Freq2,20*np.ones_like(Freq2))/20
ImgFreq=cmap(FreqN)
ImgFreq2=ImgFreq*chaosMask[:,:,None]
ImgFreq3=ImgFreq*nonchaosMask[:,:,None]

LyapN=np.minimum(Lyap2,0.006*np.ones_like(Lyap2))/0.006
ImgLyap=cmap(LyapN)
ImgLyap2=ImgLyap*chaosMask[:,:,None]
ImgLyap3=ImgLyap*nonchaosMask[:,:,None]

aux = [ImgFreq, ImgFreq2,ImgFreq3,ImgLyap,ImgLyap2,ImgLyap3, extent, ratio]


cmap=plt.get_cmap('jet')
cmap.set_under('k',1)
norm=colors.Normalize(vmin=0,vmax=0.006)



lyap_min=0
lyap_max=0.006
fig = plt.figure(1,figsize=(10,6.9))
plt.clf()
gs = gridspec.GridSpec(80, 80,wspace=0.2,hspace=0.35)
#%%subplot 3 of choosing the firing rate for chaos
ax1=plt.subplot(gs[1:24,0:19])
im1=plt.imshow(aux[1],origin='lower', interpolation='none',vmin =0, vmax = 20, extent=aux[6],aspect='auto',cmap=cmap)
#plt.xlabel(r'$\mathsf{ g_{sd} \; (mS/cm^2)}$',labelpad=0.3,fontsize='large')
plt.ylabel(r'$\mathsf{ g_{h} \; (mS/cm^2)}$',fontsize='medium')
plt.xlim([0.17,0.33])
plt.ylim([0,0.595])
plt.xticks(np.linspace(0.17,0.33,5))
plt.yticks([0,0.2,0.4,0.595],['0.0','0.2','0.4','0.6'])
ax1.text(0.177,0.61,'Selected chaotic region',fontsize='medium')
#ax1.text(0.29,0.66,'FR = %s to %s'%(minFR,maxFR),fontsize='medium')
plt.text(0.11,0.6,'A',fontsize='x-large')
#
##%%subplots 2
ax2=plt.subplot(gs[1:24,22:41])
im2=plt.imshow(aux[2],origin='lower', interpolation='none',vmin =0, vmax = 20, extent=aux[6],aspect='auto',cmap=cmap)
#plt.xlabel(r'$\mathsf{ g_{sd} \; (mS/cm^2)}$',labelpad=0.3,fontsize='large')
#    plt.ylabel(r'$\mathsf{ g_{h} \; (mS/cm^2)}$',fontsize=15)
plt.xlim([0.17,0.33])
plt.ylim([0,0.595])
plt.xticks(np.linspace(0.17,0.33,5))
plt.yticks([0,0.2,0.4,0.595],['0.0','0.2','0.4','0.6'])
ax2.text(0.175,0.61,'Selected non-chaotic region',fontsize='medium')


#%%subplot of choosing distribution for chaos condition
ax3=plt.subplot(gs[28:51,0:19])
im3=plt.imshow(aux[4],origin='lower', interpolation='none',vmin =lyap_min, vmax = lyap_max,
               extent=aux[6],aspect='auto',cmap=cmap)
plt.xlabel(r'$\mathsf{ g_{sd} \; (mS/cm^2)}$',labelpad=0.3,fontsize='medium')
plt.ylabel(r'$\mathsf{ g_{h} \; (mS/cm^2)}$',fontsize='medium')
plt.xlim([0.17,0.33])
plt.ylim([0,0.595])
plt.xticks(np.linspace(0.17,0.33,5))
plt.yticks([0,0.2,0.4,0.595],['0.0','0.2','0.4','0.6'])
#ax3.text(0.22,0.62,'chaotic',fontsize='medium')
#plt.text(0.11,0.6,'(b)',fontsize='large')
#ax1.text(0.3,0.7,'FR = 3 to 4.5',fontsize='x-large')
#    cbar=plt.colorbar(extend='max',ticks=np.linspace(0,1.6,5), pad = 0.03)
#    cbar.set_label('Lyapunov  Exponent ',fontsize=15)

#%%subplot of choosing distibution for nonchaos
ax4=plt.subplot(gs[28:51,22:41])
im4=plt.imshow(aux[5],origin='lower', interpolation='none',vmin =lyap_min, vmax = lyap_max,
           extent=aux[6],aspect='auto',cmap=cmap)
plt.xlabel(r'$\mathsf{ g_{sd} \; (mS/cm^2)}$',labelpad=0.3,fontsize='medium')
#    plt.ylabel(r'$\mathsf{ g_{h} \; (mS/cm^2)}$',fontsize=15)
plt.xlim([0.17,0.33])
plt.ylim([0,0.595])
plt.xticks(np.linspace(0.17,0.33,5))
plt.yticks([0,0.2,0.4,0.595],['0.0','0.2','0.4','0.6'])
#ax4.text(0.21,0.62,'nonchaotic',fontsize='medium')


#%%subplots 6 of the same distribution of firing rate for chaos and nonchaos neural networks
labels=['Non-chaotic','Chaotic']
ax6=plt.subplot(gs[2:44,52:80])
plt.hist((nonchaosTable2[:,-2],chaosTable[:,-2]),bins=15,color=('g','r'),label=labels)
plt.legend(fontsize='x-small',frameon=False)
plt.xlabel('Firing Rate (spikes/s)',fontsize='medium')
plt.ylabel('Event Count',fontsize='medium')
plt.xlim([minFR,maxFR])
plt.yticks(np.linspace(0.0,100,5))
plt.xticks(np.linspace(minFR,maxFR,4))
plt.text(1.8,102,'B',fontsize='x-large')
#

#setting the colorbar for the Firing rate and MLE
cax1 = fig.add_axes([0.55, 0.715, 0.01, 0.22])
cbar1=fig.colorbar(im1, extend='max',cax=cax1)
#cbar1.set_label(' Firing Rate ',fontsize='medium')
cbar1.set_label(r'FR ',fontsize='medium',labelpad=-16, y=1.15, rotation=0)
cbar1.set_ticks(np.linspace(0,20,5))
##change the appearance of ticks anf tick labbel
cbar1.ax.tick_params(labelsize='medium')


cax2 = fig.add_axes([0.55, 0.41, 0.01, 0.22])
cbar2=fig.colorbar(im3, extend='both',cax=cax2)
cbar2.set_label(r'MLE ',fontsize='medium',labelpad=-23, y=1.15, rotation=0)
cbar2.set_ticks(np.linspace(0,0.006,4))
##change the appearance of ticks anf tick labbel
cbar2.ax.tick_params(labelsize='small')


ax7=plt.subplot(gs[57:79,0:22])
plt.plot(nonchaos[:,0][1],np.mean(nonchaos[:,1:6],axis=1)[1],'g-s',label='Non-chaotic networks ')
plt.plot(nonchaos[:,0][5:],np.mean(nonchaos[:,1:6],axis=1)[5:],'g-s')

plt.plot(chaos[:,0][1],np.mean(chaos[:,1:6],axis=1)[1],'r--*',label='Chaotic networks')
plt.plot(chaos[:,0][5:],np.mean(chaos[:,1:6],axis=1)[5:],'r--*')

plt.xscale('log')
plt.xlabel(r'g', fontsize = 'medium')
plt.ylabel(r'R', fontsize='medium')

locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(0.000001,1.0 ))
ax7.xaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
ax7.xaxis.set_minor_locator(locmin)
ax7.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

y1,y2=ax7.get_ylim()
ax7.bar(3e-6,0.5,3e-6,-0.01,color='white',lw=5,clip_on=False,zorder=5)
ax7.set_ylim((y1,y2))

ax7.set_xlim((5e-7,1.5))
ax7.set_xticks((1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1))
ax7.set_xticklabels(("$0$","","$10^{-4}$","","$10^{-2}$","","$10^0$",))

t1=ax7.xaxis.get_major_ticks()[0]
t1.set_pad(5.5)

for t in ax7.xaxis.get_minor_ticks()[:27]:
    t.set_visible(False)
ax7.tick_params(labelsize=9)





ax8=plt.subplot(gs[57:79,29:51])
plt.plot(nonchaos[:,0][1],np.mean(nonchaos[:,11:16],axis=1)[1],'g-s')
plt.plot(nonchaos[:,0][5:],np.mean(nonchaos[:,11:16],axis=1)[5:],'g-s')

plt.plot(chaos[:,0][1],np.mean(chaos[:,11:16],axis=1)[1],'r--*')
plt.plot(chaos[:,0][5:],np.mean(chaos[:,11:16],axis=1)[5:],'r--*')

plt.xscale('log')
plt.xlabel(r'g', fontsize = 'medium')
plt.ylabel(u'Metastability', fontsize='medium')
plt.ylim([-0.0015,0.031])
plt.yticks(np.linspace(0,0.03,4))

locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(0.000001,1.0 ))
ax8.xaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
ax8.xaxis.set_minor_locator(locmin)
ax8.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

y1,y2=ax8.get_ylim()
ax8.bar(3e-6,0.11,3e-6,-0.1,color='white',lw=5,clip_on=False,zorder=5)
ax8.set_ylim((y1,y2))

ax8.set_xlim((5e-7,1.5))
ax8.set_xticks((1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1))
ax8.set_xticklabels(("$0$","","$10^{-4}$","","$10^{-2}$","","$10^0$",))

t1=ax8.xaxis.get_major_ticks()[0]
t1.set_pad(5.5)

for t in ax8.xaxis.get_minor_ticks()[:27]:
    t.set_visible(False)
ax8.tick_params(labelsize=9)


ax9=plt.subplot(gs[57:79,58:80])
plt.plot(nonchaos[:,0][1],np.mean(nonchaos[:,6:11],axis=1)[1],'g-s')
plt.plot(nonchaos[:,0][5:],np.mean(nonchaos[:,6:11],axis=1)[5:],'g-s')

plt.plot(chaos[:,0][1],np.mean(chaos[:,6:11],axis=1)[1],'r--*')
plt.plot(chaos[:,0][5:],np.mean(chaos[:,6:11],axis=1)[5:],'r--*')

plt.xscale('log')
plt.ylim([-0.001,0.008])
plt.xlabel(r'g', fontsize = 'medium')
plt.ylabel(r'MLE',fontsize='medium')
plt.yticks([0,0.002,0.004,0.006,0.008],['0.000','0.002','0.004','0.006','0.008'])

locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(0.000001,1.0 ))
ax9.xaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )) 
ax9.xaxis.set_minor_locator(locmin)
ax9.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

y1,y2=ax9.get_ylim()
ax9.bar(3e-6,0.1,3e-6,-0.1,color='white',lw=5,clip_on=False,zorder=5)
ax9.set_ylim((y1,y2))

ax9.set_xlim((5e-7,1.5))
ax9.set_xticks((1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1))
ax9.set_xticklabels(("$0$","","$10^{-4}$","","$10^{-2}$","","$10^0$",))

t1=ax9.xaxis.get_major_ticks()[0]
t1.set_pad(5.5)

for t in ax9.xaxis.get_minor_ticks()[:27]:
    t.set_visible(False)
ax9.tick_params(labelsize=9)





#
fig.subplots_adjust(left=0.1,bottom=0.08,top=0.98,right=0.97)
plt.savefig('S1_Fig_ghgsdDis.png',dpi=300)
#S1_Fig_ghgsdDis
#
#
#
