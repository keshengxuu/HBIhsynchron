# ranges sets 
#range1
ghmin1, ghmax1, gsdmin1, gsdmax1= 0.35, 0.4, 0.197, 0.207 #range chaos 1 gh/gsd T36
ghmin2, ghmax2, gsdmin2, gsdmax2= 0.1, 0.15, 0.197, 0.207 #range nonchaos 1 gh/gsd T36

#range2
ghmin1b, ghmax1b, gsdmin1b, gsdmax1b= 0.35, 0.40, 0.275, 0.285 #range chaos 4 gh/gsd T36
ghmin2b, ghmax2b, gsdmin2b, gsdmax2b= 0.18, 0.23, 0.275,0.285 #range nonchaos 4 gh/gsd T36


# paremeters for the network
Nnode = 50     # the number of neurons
lyap_min=0         #minmum lyapunov value
lyap_max=1.6         #maxsimum lyapunov value
dt = 0.05       # the time step of simulation
Xdata=0
Ydata=4
Ydata2=3
tEnd=50000
nsim=50
tbin=40
rseed0=[0,20,25,33,40]  # the seed for simulation
ranges0=[2,2,2,2,2] # the index of ranges2
Nindex=1    # The number of  example ranges for fixed size
