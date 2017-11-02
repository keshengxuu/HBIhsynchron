# -*- coding: utf-8 -*-
"""
Firing rate statistics
FIXED DISTRIBUTION OF gsd, gh

"""
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
os.chdir(parentdir)

from lib import loaddata as ld
from lib import Figures_plot as Fig
from lib import filter_function as fc

#%%
#Table=np.loadtxt("spikes_gsdgsrT36gh06/finalTable.txt",skiprows=1)
execfile('lib/Defaultsets.py')

#%%
# calculating the mean fraction of  bursting events
sp_chao,sp_nonchao=ld.load_data_spike(rseed0,ranges0)



#for i in range(5):
Burst_chaos, Burst_nonch = fc.frac_burst(sp_chao,sp_nonchao, Nindex)
#


# call the figs for plotting
# elem and elemb are  the  bursting events for chaotic and  nonchaotic,respectively.
Fig.FigX_plots(Burst_chaos, Burst_nonch)


