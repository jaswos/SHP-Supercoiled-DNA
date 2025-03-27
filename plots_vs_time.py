

import numpy as np
import matplotlib.pyplot as plt#
import sys


fileLen= 602  # length of file in lines. First line is not read
cm = 1/2.54  # centimeters in inches
Tmax= 30e6    # duration of simulation in steps
porsion_list= [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
sample_interval = 50000
cutoff = int(4e6/sample_interval)

run=3  # only look at one run

########################################################
## load some data sets

data = {}  # store data in a dict

output_name = sys.argv[1]

for porsion in porsion_list:
    data[ (porsion,run) ]=np.zeros((fileLen-1,10))
    rf = open(f'split_measurements/scMeasure_L300_porsion{porsion}_run{run}.dat', 'r')
    for i in range(fileLen):
        line = rf.readline()
        if i>0:
            line = line.split()
            for j in range(len(line)):
                data[(porsion,run)][i-1][j]=np.float(line[j])
    print(f"Min/Max Stiff Writhe: {np.min(data[(porsion,run)][:,8])}, {np.max(data[(porsion,run)][:,8])}")
    print(f"Min/Max Normal Writhe: {np.min(data[(porsion,run)][:,9])}, {np.max(data[(porsion,run)][:,9])}")

    rf.close()

# so data[(porsion,run)] is a numpy array with 10 columns
#                      data[(porsion,run)][:,0]   are times in steps
#                      data[(porsion,run)][:,1]   is normal radius of gyration
#                      data[(porsion,run)][:,2]   is stiff radius of gyration
#                      data[(porsion,run)][:,3]   is total twist
#                      data[(porsion,run)][:,4]   is normal twist
#                      data[(porsion,run)][:,5]   is stiff twist
#                      data[(porsion,run)][:,6]   is twist deficit
#                      data[(porsion,run)][:,7]   is total writhe
#                      data[(porsion,run)][:,8]   is stiff writhe
#                      data[(porsion,run)][:,9]   is normal writhe


########################################################
## Now do some plots as a function of time

plt.rcParams.update({'font.size': 6})

fig, axs = plt.subplots(1, 5)
fig.set_figwidth(60*cm)
fig.set_figheight(12*cm)

for porsion in porsion_list:
    axs[0].plot(data[(porsion,run)][cutoff:,0]/1e6,data[(porsion,run)][cutoff:,1], label=f'porsion={porsion}', 
                linewidth = 0.8)
axs[0].set_xlim(cutoff * sample_interval/1e6, Tmax/1e6)
axs[0].set_xlabel('time (10^6 steps)')
axs[0].set_ylabel('radius of gyration')
axs[0].legend(loc='upper right')

for porsion in porsion_list:
    axs[1].plot(data[(porsion,run)][cutoff:,0]/1e6,data[(porsion,run)][cutoff:,3], label=f'porsion={porsion}',
                linewidth = 0.4)
axs[1].set_xlim(cutoff * sample_interval/1e6, Tmax/1e6)
axs[1].set_xlabel('time (10^6 steps)')
axs[1].set_ylabel('Total Twist')
axs[1].legend(loc='upper right')

for porsion in porsion_list:
    axs[2].plot(data[(porsion,run)][cutoff:,0]/1e6,data[(porsion,run)][cutoff:,5], label=f'porsion={porsion}',
                linewidth = 0.4)
axs[2].set_xlim(cutoff * sample_interval/1e6, Tmax/1e6)
axs[2].set_xlabel('time (10^6 steps)')
axs[2].set_ylabel('Stiff Twist')
axs[2].legend(loc='upper right')

for porsion in porsion_list:
    axs[3].plot(data[(porsion,run)][cutoff:,0]/1e6,data[(porsion,run)][cutoff:,8], label=f'porsion={porsion}',
                linewidth = 0.4)
axs[3].set_xlim(cutoff * sample_interval/1e6, Tmax/1e6)
axs[3].set_xlabel('time (10^6 steps)')
axs[3].set_ylabel('Stiff Writhe')
axs[3].legend(loc='upper right')

for porsion in porsion_list:
    axs[4].plot(data[(porsion,run)][cutoff:,0]/1e6,data[(porsion,run)][cutoff:,9], label=f'porsion={porsion}',
                linewidth = 0.4)
axs[4].set_xlim(cutoff * sample_interval/1e6, Tmax/1e6)
axs[4].set_xlabel('time (10^6 steps)')
axs[4].set_ylabel('Normal Writhe')
axs[4].legend(loc='upper right')



plt.savefig(f"plots/{output_name}.pdf", format="pdf", bbox_inches="tight", dpi = 400)