

import numpy as np
import matplotlib.pyplot as plt
# the lesser known counterpart to skegness
import scipy.stats
import sys
import seaborn as sns

cm = 1/2.54  # centimeters in inches
Tmax=14e6    # duration of simulation in steps
tlist= np.array([1, 2, 3, 4, 5, 10])
runlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

def moving_average(array, window_size):
    return np.convolve(array, np.ones(window_size)/window_size, mode='valid')



if not len(sys.argv) == 2:
    print("Command usage: python plots_decay.py which_lk")
    print("where   which_lk         is an string, initial linking number of polymer")


which_lk = sys.argv[1]
########################################################
## load some data sets

data = {}  # store data in a dict

checkLines = open(f"measurements/scMeasure_L300_t{tlist[0]}_run{runlist[0]}.dat")
fileLen = len(checkLines.readlines())
checkLines.close()


for t in tlist:
    for run in runlist:
        data[(t, run)] = np.zeros((fileLen-1,5))
        rf = open(f"measurements/scMeasure_L300_t{t}_run{run}.dat", 'r')
        for i in range(fileLen):
            line = rf.readline()
            if i>0:
                line = line.split()
                for j in range(len(line)):
                    data[(t, run)][i-1][j]=float(line[j])
        rf.close()

# so data[(t,run)] is a numpy array with five columns
#                      data[(t,run)][:,0]   are times in steps
#                      data[(t,run)][:,1]   is radius of gyration
#                      data[(t,run)][:,2]   is twist
#                      data[(t,run)][:,3]   is writhe 
#                      data[(t,run)][:,5]   is linking number


## Decay anaylsis; this is a non-equilibrium process, so mainly
## looking at the decay timescales

torque_index = int(2.2e6/10000)
decay_end = int(4.9e6/10000)
mean_from = int(6e6/10000)

decorrelationTime = 57250

mean_twists, std_twists = [], []

for t in tlist:
    twist_means = []
    std_means = []
    for run in runlist:
        twist_means.append(np.mean(data[(t, run)][mean_from:, 2]))
        std_means.append(np.std(data[(t, run)][mean_from:, 2]))
    mean_twists.append(np.mean(twist_means))
    std_twists.append(np.mean(std_means)/(np.sqrt((Tmax - (10000*mean_from))/decorrelationTime)))

print(f"Total samples: {(Tmax - (10000*mean_from))/decorrelationTime}")
    
lk_decay = np.zeros((len(tlist), decay_end - torque_index))

for i in range(len(tlist)):
    for j in range(len(lk_decay[0,:])):
        mean_lk_at_timestep = 0
        for k in range(len(runlist)):
            #print(f"doing t {tlist[i]}, timestep {j}, run {runlist[k]}")
            mean_lk_at_timestep += data[(tlist[i], runlist[k])][torque_index + j, 4]
        lk_decay[i, j] = mean_lk_at_timestep/(len(runlist))
        
all_quasi_equilibrium_twist = []
for i in range(len(tlist)):
    all_quasi_equilibrium_twist.append([])

for i in range(len(tlist)):
    for j in range(len(runlist)):
        for k in range(399):
            all_quasi_equilibrium_twist[i].append(data[tlist[i], runlist[j]][int(1000 + k), 2])
        
lk_decay_slopes = []
for i in range(len(tlist)):
    lk_decay_slopes.append(scipy.stats.linregress(np.log(data[(tlist[i], runlist[0])][torque_index:decay_end, 0]), np.log(lk_decay[i, :])).slope)

#Trying to plot the mean slop of the log log decay

times = data[(tlist[0], runlist[0])][torque_index:decay_end, 0]
log_times = np.log(times)

mean_slope = np.mean(lk_decay_slopes)

mean_start_value = 0
for i in range(len(tlist)):
    mean_start_value += lk_decay[i, 0]
mean_start_value /= len(tlist)
logged_mean_start = np.log(mean_start_value)

mean_decay = [logged_mean_start]

for i in range(len(log_times)-1):
    mean_decay.append(mean_decay[-1] + (mean_slope)*(log_times[i+1] - log_times[i]))
    
colour = ['tab:red', 'tab:orange', 'tab:green', 'tab:cyan', 'tab:blue', 'tab:purple']

plt.rcParams.update({'font.size': 14})

"""
fig, axs = plt.subplots(1, 6)
fig.set_figwidth(110*cm)
fig.set_figheight(10*cm)

axs[0].plot(tlist, mean_twists)
axs[0].errorbar(tlist, mean_twists, yerr=std_twists, fmt='o', capsize=5)
axs[0].set_xlabel('Torque applied to polymer / k_B T')
axs[0].set_ylabel('Mean Twist across all runs')

axs[1].plot(tlist, std_twists)
axs[1].set_xlabel('Torque applied to polymer / k_B T')
axs[1].set_ylabel('Standard deviation on the twist in the nicked regime')

for i in range(len(tlist)):
    axs[2].plot(log_times, np.log(lk_decay[i, :]), label = f't={tlist[i]}', alpha = 0.8, color = colour[i])
axs[2].plot(log_times, mean_decay, label = "Mean slope", linewidth = 2, linestyle = 'dashed', color = 'k')
axs[2].set_xlabel('Log Timestep')
axs[2].set_ylabel('Log of averaged writhe across all runs')
axs[2].legend(loc = 'lower left')

for i in range(len(tlist)):
    axs[3].plot(times, lk_decay[i, :], label = f't={tlist[i]}', color = colour[i])
axs[3].set_xlabel('Timestep')
axs[3].set_ylabel('Averaged writhe across all runs')
axs[3].legend(loc = 'upper right')

print(f"Average decay slope: {np.mean(lk_decay_slopes)}")
print(f"Std of the decay: {np.std(mean_twists)}")

#Histogramming
N = 5
def cleaner(array):
    cleaned = []
    for i in range(len(array)):
        val = abs(array[i])
        if val > 0.00001:
            cleaned.append(array[i])
    return cleaned


#Plot 2 different series, using indices to select

plot1_index = 0
plot2_index = 5

sns.kdeplot(cleaner(all_quasi_equilibrium_twist[plot1_index]), bw_adjust=0.4, color=colour[plot1_index], linewidth=2, ax=axs[4], label = f"torque: {tlist[plot1_index]}")
sns.kdeplot(cleaner(all_quasi_equilibrium_twist[plot2_index]), bw_adjust=0.4, color=colour[plot2_index], linewidth=2, ax=axs[4], label = f"torque: {tlist[plot2_index]}")
sns.histplot(cleaner(all_quasi_equilibrium_twist[0]), binwidth=0.05, kde=True, alpha = 0.5, color=colour[0], ax=axs[4])
sns.histplot(cleaner(all_quasi_equilibrium_twist[5]), binwidth=0.05, kde=True, alpha = 0.5, color=colour[5], ax=axs[4])

axs[4].axvline(0, color='lightgray', linestyle='dotted', linewidth=2.0)
axs[4].set_title(f"Twist distribution over the quasi-equilibrium (>{mean_from/100}e6) regime")
axs[4].set_xlabel("Twist Value")
axs[4].set_ylabel("Frequency")
axs[4].legend()


for i in range(len(tlist)):
    sns.kdeplot(cleaner(all_quasi_equilibrium_twist[i]), linewidth=1, ax=axs[5], label = f'torque: {tlist[i]}', color = colour[i])
axs[5].set_title("KDE fits to histograms of twist")
axs[5].axvline(0, color='lightgray', linestyle='dotted', linewidth=2.0)
axs[5].set_xlabel("Twist Value")
axs[5].set_ylabel("Density")
axs[5].legend(loc = 'upper right')

plt.savefig(f"plots/nonEquilibrium_lk{which_lk}_all.pdf", dpi = 600, format="pdf", bbox_inches="tight")


#Producing some desired plots

slope, intercept = np.polyfit(tlist[:-1], mean_twists[:-1], 1)  # 1st-degree polynomial (linear fit)
mean_twist_fit = slope * tlist + intercept  # Line equation: y = mx + c


mean_writhe, std_writhe = [], []

for t in tlist:
    writhe_means = []
    std_means = []
    for run in runlist:
        writhe_means.append(np.mean(data[(t, run)][mean_from:, 3]))  # Column 3 = writhe
        std_means.append(np.std(data[(t, run)][mean_from:, 3]))
    mean_writhe.append(np.mean(writhe_means))
    std_writhe.append(np.mean(std_means) / np.sqrt((Tmax - (10000 * mean_from)) / decorrelationTime))

# Perform linear regression for writhe
slope_writhe, intercept_writhe = np.polyfit(tlist[:-1], mean_writhe[:-1], 1)
mean_writhe_fit = slope_writhe * tlist + intercept_writhe



# Plot twist and writhe with regression lines
fig, axs = plt.subplots(1, 3)
fig.set_figwidth(34 * cm)
fig.set_figheight(10 * cm)

# Plot Writhe
axs[0].errorbar(tlist, mean_writhe, yerr=std_writhe, fmt='s', capsize=5, color='tab:red', label="Mean Writhe")  # Data points + error bars
axs[0].plot(tlist[:-1], mean_writhe_fit[:-1], color='tab:orange', label="Writhe Fit")  # Regression line
axs[0].plot(tlist[-2:], mean_writhe_fit[-2:], color='tab:orange', linestyle='dashed')  # Dashed extrapolated line


axs[0].errorbar(tlist, mean_twists, yerr=std_twists, fmt='o', capsize=5, color='tab:cyan', label="Mean Twist")  # Data points + error bars
axs[0].plot(tlist[:-1], mean_twist_fit[:-1], color='tab:blue', label="Twist Fit")  # Regression line
axs[0].plot(tlist[-2:], mean_twist_fit[-2:], color='tab:blue', linestyle='dashed')  # Dashed extrapolated line


# Custom Legend with Markers and Lines
from matplotlib.lines import Line2D

legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:cyan', markersize=8, label='Mean Twist'),
    Line2D([0], [0], linestyle='-', color='tab:blue', lw=2, label='Twist Fit'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor='tab:red', markersize=8, label='Mean Writhe'),
    Line2D([0], [0], linestyle='-', color='tab:orange', lw=2, label='Writhe Fit'),
]

axs[0].legend(handles=legend_elements)

axs[0].set_xlabel(r'Torque applied to polymer  ($k_BT$)')
axs[0].set_ylabel('Mean Twist & Writhe across all runs')
axs[0].legend()

print(mean_twists)
print(std_twists)


sns.kdeplot(cleaner(all_quasi_equilibrium_twist[plot1_index]), bw_adjust=0.4, color=colour[plot1_index], linewidth=2, ax=axs[1], label = f"t={tlist[plot1_index]}")
sns.kdeplot(cleaner(all_quasi_equilibrium_twist[plot2_index]), bw_adjust=0.4, color=colour[plot2_index], linewidth=2, ax=axs[1], label = f"t={tlist[plot2_index]}")
sns.histplot(cleaner(all_quasi_equilibrium_twist[0]), binwidth=0.05, kde=True, alpha = 0.5, color=colour[0], ax=axs[1])
sns.histplot(cleaner(all_quasi_equilibrium_twist[5]), binwidth=0.05, kde=True, alpha = 0.5, color=colour[5], ax=axs[1])

axs[1].axvline(0, color='lightgray', linestyle='dotted', linewidth=2.0)
axs[1].set_title(f"Twist distribution")
axs[1].set_xlabel("Twist Value")
axs[1].set_ylabel("Frequency")
axs[1].legend()


for i in range(len(tlist)):
    sns.kdeplot(cleaner(all_quasi_equilibrium_twist[i]), linewidth=1, ax=axs[2], label = f't={tlist[i]}', color = colour[i])
axs[2].set_title("KDE fits to histograms of twist")
axs[2].axvline(0, color='lightgray', linestyle='dotted', linewidth=2.0)
axs[2].set_xlabel("Twist Value")
axs[2].set_ylabel("Density")
axs[2].legend(loc = 'upper right')

plt.savefig(f"plots/twist_vs_torque_{which_lk}.png", dpi = 600, format="png", bbox_inches="tight")
"""

fig, axs = plt.subplots()
fig.set_figwidth(15*cm)
fig.set_figheight(12*cm)

smooth = 5

for i in range(len(tlist)):
    axs.plot(log_times[smooth-1:], moving_average(np.log(lk_decay[i, :]), smooth), label = f'T={tlist[i]}', alpha = 0.8, color = colour[i])
axs.plot(log_times, mean_decay, label = "Mean slope", linewidth = 2, linestyle = 'dashed', color = 'k')
axs.set_xlabel('Log Timestep')
axs.set_ylabel(r'ln$\left(\overline{Lk}\right)$')
axs.legend(loc = 'lower left')

plt.savefig(f"plots/for_report.png", dpi = 600, format="png", bbox_inches="tight")
