# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 16:28:46 2025

@author: hp
"""

import numpy as np
import matplotlib.pyplot as plt


fileLen= 402  # length of file in lines. First line is not read
cm = 1/2.54  # centimeters in inches
Tmax= 20e6    # duration of simulation in steps
bs_list= [20.0, 22.0, 24.0, 26.0, 28.0, 30.0]

# The data sampling follows the dumpfile, which has data every 50000 timesteps
# As the data is assumed to be unphysical before 4 million timesteps, this 
# index allows us to only plot physically significant data
cutoff_index = int(4000000/50000)

run=1  # only look at one run

########################################################
## load some data sets

data = {}  # store data in a dict

for bs in bs_list:
    data[ (bs,run) ]=np.zeros((fileLen-1,5))
    rf = open(f'measurements/scMeasure_L300_bs{bs}_run{run}.dat', 'r')
    for i in range(fileLen):
        line = rf.readline()
        if i>0:
            line = line.split()
            for j in range(len(line)):
                data[(bs,run)][i-1][j]=np.float(line[j])
    rf.close()

# so data[(bs,run)] is a numpy array with five columns
#                      data[(bs,run)][:,0]   are times in steps
#                      data[(bs,run)][:,1]   is radius of gyration
#                      data[(bs,run)][:,2]   is twist
#                      data[(bs,run)][:,3]   is writhe 
#                      data[(bs,run)][:,5]   is linking number


########################################################
## Now do some plots as a function of time

plt.rcParams.update({'font.size': 6})

fig, axs = plt.subplots(1, 4)
fig.set_figwidth(45*cm)
fig.set_figheight(12*cm)

for bs in bs_list:
    axs[0].plot(data[(bs,run)][cutoff_index:,0]/1e6,data[(bs,run)][cutoff_index:,1], label=f'bs={bs}', 
                linewidth = 0.8)
axs[0].set_xlim(0, Tmax/1e6)
axs[0].set_xlabel('time (10^6 steps)')
axs[0].set_ylabel('radius of gyration')
axs[0].legend(loc='upper right')

for bs in bs_list:
    axs[1].plot(data[(bs,run)][cutoff_index:,0]/1e6,data[(bs,run)][cutoff_index:,2], label=f'bs={bs}',
                linewidth = 0.8)
axs[1].set_xlim(0, Tmax/1e6)
axs[1].set_xlabel('time (10^6 steps)')
axs[1].set_ylabel('twist')
axs[1].legend(loc='upper right')

for bs in bs_list:
    axs[2].plot(data[(bs,run)][cutoff_index:,0]/1e6,data[(bs,run)][cutoff_index:,3], linewidth = 0.8)
axs[2].set_xlim(0, Tmax/1e6)
axs[2].set_xlabel('time (10^6 steps)')
axs[2].set_ylabel('writhe')

plt.savefig("plots/vsTime_run1_higher.pdf", format="pdf", bbox_inches="tight", dpi = 400)
