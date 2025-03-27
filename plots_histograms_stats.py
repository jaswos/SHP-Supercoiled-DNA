

import numpy as np
import matplotlib.pyplot as plt

fileLen=402  # length of file in lines. First line is not read
cm = 1/2.54  # centimeters in inches
Tmax=30e6    # duration of simulation in steps
cutFirst=4e6  # ignore the first this many steps, as it has not reached a steady state 

bslist=[10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0]
# 
runs=list(range(1, 11))


####################################
## load data sets

data = {}  # store data in a dict

for bs in bslist:
    for run in runs:
        data[ (bs,run) ]=np.zeros((fileLen-1,10))
        rf = open(f'split_measurements/scMeasure_L300_bs{bs}_run{run}.dat', 'r')
        for i in range(fileLen):
            line = rf.readline()
            if i>0:
                line = line.split()
                for j in range(len(line)):
                    data[(bs,run)][i-1][j]=float(line[j])
        rf.close()

# so data[(bs,run)] is a numpy array with 10 columns
#                      data[(bs,run)][:,0]   are times in steps
#                      data[(bs,run)][:,1]   is normal radius of gyration
#                      data[(bs,run)][:,2]   is stiff radius of gyration
#                      data[(bs,run)][:,3]   is total twist
#                      data[(bs,run)][:,4]   is normal twist
#                      data[(bs,run)][:,5]   is stiff twist
#                      data[(bs,run)][:,6]   is twist deficit
#                      data[(bs,run)][:,7]   is total writhe
#                      data[(bs,run)][:,8]   is stiff writhe
#                      data[(bs,run)][:,9]   is normal writhe
####################################


## collect all data, removing first section from each sim

allNRg = {}  # a dictionary which will contain a big np array for each bs
allSRg = {}
allNTw = {}
allSTw = {}
allNW = {}
allSW = {}

removeTo=int(4000000/50000)   # find the index where the time = cutFirst

"""
this has been entered manually for now because it's half 1, i want to sleep,
and i don't knwo how many more errors i can handle before i start crying...
"""

for bs in bslist:
    # run 1
    allNRg[bs] = data[(bs,1)][removeTo:,1]
    allSRg[bs] = data[(bs,1)][removeTo:,2]
    allNTw[bs] = data[(bs,1)][removeTo:,3]
    allSTw[bs] = data[(bs,1)][removeTo:,4]
    allNW[bs] = data[(bs,1)][removeTo:,9]  # Normal Writhe
    allSW[bs] = data[(bs,1)][removeTo:,8]  # Stiff Writhe
    
    # subsequent runs - append them onto this array
    for run in runs[1:]:
        allNRg[bs] = np.append(allNRg[bs], data[(bs,run)][removeTo:,1])
        allSRg[bs] = np.append(allSRg[bs], data[(bs,run)][removeTo:,2])
        allNTw[bs] = np.append(allNTw[bs], data[(bs,run)][removeTo:,3])
        allSTw[bs] = np.append(allSTw[bs], data[(bs,run)][removeTo:,4])
        allNW[bs] = np.append(allNW[bs], data[(bs,run)][removeTo:,9])
        allSW[bs] = np.append(allSW[bs], data[(bs,run)][removeTo:,8])

# Define bin widths
NRg_bw, SRg_bw = 1.20, 1.20
NTw_bw, STw_bw = 0.05, 0.05  # Twist fraction
NW_bw, SW_bw = 0.2, 0.2   # Writhe fraction

# Initialize dictionaries
binnedNRg, binnedSRg, binnedNTw, binnedSTw = {}, {}, {}, {}
binnedNW, binnedSW = {}, {}  # New dictionaries for writhe histograms

for bs in bslist:
    for dataset, bw, binned_dict in [(allNRg, NRg_bw, binnedNRg), (allSRg, SRg_bw, binnedSRg),
                                     (allNTw, NTw_bw, binnedNTw), (allSTw, STw_bw, binnedSTw),
                                     (allNW, NW_bw, binnedNW), (allSW, SW_bw, binnedSW)]:
        bins = {}
        for value in dataset[bs]:
            bins[np.round(value / bw) * bw] = bins.get(np.round(value / bw) * bw, 0) + 1
        binned_dict[bs] = np.array(list(bins.items()))
        binned_dict[bs] = binned_dict[bs][binned_dict[bs][:, 0].argsort()]

plt.rcParams.update({'font.size': 6})
fig, axs = plt.subplots(1, 6)
fig.set_figwidth(57 * cm)
fig.set_figheight(10 * cm)

plot_info = [
    (binnedNRg, 'radius of gyration normal section'),
    (binnedSRg, 'radius of gyration stiff section'),
    (binnedNTw, 'twist normal section'),
    (binnedSTw, 'twist stiff section'),
    (binnedNW, 'normal writhe'),
    (binnedSW, 'stiff writhe')
]

for i, (binned_dict, xlabel) in enumerate(plot_info):
    for bs in bslist:
        axs[i].plot(binned_dict[bs][:, 0], binned_dict[bs][:, 1], marker='o', markersize=2,
                    label=f'bs={bs}', linewidth=0.8)
    axs[i].set_xlabel(xlabel)
    axs[i].set_ylabel('counts')
    which_to_show = [0, 1, 2, 3]
    if i in which_to_show:
        axs[i].legend(loc='upper right')

plt.savefig("histograms_extended.pdf", format="pdf", bbox_inches="tight")


##################################
# Calculate some statistics 

# Initialize dictionaries to store mean and standard deviation
import numpy as np
import matplotlib.pyplot as plt

# Dictionaries to store means and standard deviations
meanNRg, stdNRg = {}, {}
meanSRg, stdSRg = {}, {}
meanNTw, stdNTw = {}, {}
meanSTw, stdSTw = {}, {}
meanNW, stdNW = {}, {}  # Normal Writhe
meanSW, stdSW = {}, {}  # Stiff Writhe

for bs in bslist:
    meanNRg[bs] = np.mean(allNRg[bs])
    stdNRg[bs] = np.std(allNRg[bs])

    meanSRg[bs] = np.mean(allSRg[bs])
    stdSRg[bs] = np.std(allSRg[bs])

    meanNTw[bs] = np.mean(allNTw[bs])
    stdNTw[bs] = np.std(allNTw[bs])

    meanSTw[bs] = np.mean(allSTw[bs])
    stdSTw[bs] = np.std(allSTw[bs])

    meanNW[bs] = np.mean(allNW[bs])  # Normal Writhe
    stdNW[bs] = np.std(allNW[bs])

    meanSW[bs] = np.mean(allSW[bs])  # Stiff Writhe
    stdSW[bs] = np.std(allSW[bs])

# Convert dictionaries to arrays for plotting
bs_array = np.array(bslist)

# Compute standard error of the mean (SEM)
SEM_factor = np.sqrt(10)  # since SEM = std / sqrt(N), with N=10
NRg_sem = np.array([stdNRg[bs] / SEM_factor for bs in bslist])
SRg_sem = np.array([stdSRg[bs] / SEM_factor for bs in bslist])
NTw_sem = np.array([stdNTw[bs] / SEM_factor for bs in bslist])
STw_sem = np.array([stdSTw[bs] / SEM_factor for bs in bslist])
NW_sem = np.array([stdNW[bs] / SEM_factor for bs in bslist])
SW_sem = np.array([stdSW[bs] / SEM_factor for bs in bslist])

plt.rcParams.update({'font.size': 6})
fig, axs = plt.subplots(1, 6)
fig.set_figwidth(57 * cm)
fig.set_figheight(10 * cm)

plot_info = [
    ([meanNRg, stdNRg], 'radius of gyration normal section', 'NRg'),
    ([meanSRg, stdSRg], 'radius of gyration stiff section', 'SRg'),
    ([meanNTw, stdNTw], 'twist stiff section', 'NTw'),
    ([meanSTw, stdSTw], 'twist normal section', 'STw'),
    ([meanNW, stdNW], 'normal writhe', 'NW'),
    ([meanSW, stdSW], 'stiff writhe', 'SW')
]

SEM_factor = np.sqrt(10)  # since SEM = std / sqrt(N), with N=10

for i, ((mean_dict, std_dict), xlabel, label) in enumerate(plot_info):
    mean_values = np.array([mean_dict[bs] for bs in bslist])
    sem_values = np.array([std_dict[bs] / SEM_factor for bs in bslist])

    axs[i].errorbar(bslist, mean_values, yerr=sem_values, fmt='o-', label=label)
    axs[i].set_xlabel(xlabel)
    axs[i].set_ylabel('mean ± SEM')
    axs[i].legend(loc='upper right')

plt.savefig("means_with_errorbars.pdf", format="pdf", bbox_inches="tight")





##################################
# Plot the stats

## Will want to plot mean Tw vs bs etc.
