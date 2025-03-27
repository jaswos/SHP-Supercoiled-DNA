

import numpy as np
import matplotlib.pyplot as plt
import sys

cm = 1/2.54  # centimeters in inches
Tmax=14e6    # duration of simulation in steps
tlist= [1, 2, 3, 4, 5, 10]

run = sys.argv[1] # only look at one run
which_lk = sys.argv[2]

if not len(sys.argv) == 3:
    print("Command usage: python plots_vs_time.py run which_lk")
    print("where   run              is an string, which run files to look at")
    print("where   which_lk         is an string, initial linking number of polymer")


########################################################
## load some data sets

data = {}  # store data in a dict

checkLines = open(f"measurements/scMeasure_L300_t1_run1.dat")
fileLen = len(checkLines.readlines())
checkLines.close()


for t in tlist:
    data[ (t,run) ]=np.zeros((fileLen-1,5))
    rf = open(f"measurements/scMeasure_L300_t{t}_run{run}.dat", 'r')
    for i in range(fileLen):
        line = rf.readline()
        if i>0:
            line = line.split()
            for j in range(len(line)):
                data[(t,run)][i-1][j]=float(line[j])
    rf.close()

# so data[(t,run)] is a numpy array with five columns
#                      data[(t,run)][:,0]   are times in steps
#                      data[(t,run)][:,1]   is radius of gyration
#                      data[(t,run)][:,2]   is twist
#                      data[(t,run)][:,3]   is writhe 
#                      data[(t,run)][:,5]   is linking number


########################################################
## Now do some plots as a function of time

plt.rcParams.update({'font.size': 11})

fig, axs = plt.subplots(1, 3)
fig.set_figwidth(30*cm)
fig.set_figheight(10*cm)

for t in tlist:
    axs[0].plot(data[(t,run)][:,0]/1e6, data[(t,run)][:,2], label=f't={t}', linewidth=0.3)
axs[0].set_xlim(0, Tmax/1e6)
axs[0].set_xlabel(r'Timesteps ($10^6$)')
axs[0].set_ylabel('Twist')
axs[0].legend(loc='upper right')

for t in tlist:
    axs[1].plot(data[(t,run)][:,0]/1e6, data[(t,run)][:,3], label=f't={t}', linewidth=0.3)
axs[1].set_xlim(0, Tmax/1e6)
axs[1].set_xlabel(r'Timesteps ($10^6$)')
axs[1].set_ylabel('Writhe')
axs[1].legend(loc='upper right')

for t in tlist:
    axs[2].plot(data[(t,run)][:,0]/1e6, data[(t,run)][:,4], label=f't={t}', linewidth=0.3)
axs[2].set_xlim(0, Tmax/1e6)
axs[2].set_xlabel(r'Timesteps ($10^6$)')
axs[2].set_ylabel('Linking number')
axs[2].legend(loc='upper right')

for legend in [axs[j].legend() for j in [0, 1, 2]]:
    for line in legend.get_lines():
        line.set_linewidth(1.2)  # Set thicker legend lines

plt.savefig(f"plots/vsTime_lk{which_lk}_run{run}.png", format="png", dpi = 600, bbox_inches="tight")


plt.show()