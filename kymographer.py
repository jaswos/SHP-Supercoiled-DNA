# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 21:38:09 2025

@author: hp
"""

import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np
import sys 

if not len(sys.argv) == 7:
    print("Command usage: python ./kymographer.py Npartitions dumpInterval smudge infile colour_style plot_label")
    print("where   Npartitions       is an integer, how many partitions the polymer is separated into")
    print("        dumpInterval      is an integer, the frequency with which data was output to dumpfile")
    print("        smudge            is an integer, giving the width of the averagin bracket")
    print("        infile            is a file, processed partition data")
    print("        colour_style      is a string, check the colormap documentation for options")
    print("        plot_label        is a string, added to output file names")
    exit()

Npartitions = int(sys.argv[1])
dumpInterval = int(sys.argv[2])
smudge = int(sys.argv[3])
infile = sys.argv[4]
colour_style = str(sys.argv[5])
plot_label = sys.argv[6]

rawfile = open(infile, "r")
lines = rawfile.readlines()
Nlines = len(lines)

data = np.zeros((Nlines - 1, Npartitions))

for i in range(1, Nlines):
    # note the first two entries are the timestep and the total twist
    split_line = lines[i].split()
    for j in range(2, Npartitions+2):
        #print(f"i: {i}, j: {j}")
        data[i-1, j-2] = float(split_line[j])

data_smudge = np.zeros((Nlines - 1, Npartitions))

# smudge is a cute variable to define how wide the averaging bracket is
for i in range(0, len(data[:, 0])):
     for j in range(0, len(data[0, :])):
         to_be_averaged = []
         for k in range(-smudge, smudge):
            to_be_averaged.append(data[i, ((j + k)%Npartitions - 1)])
         data_smudge[i, j] = np.mean(to_be_averaged)
     
plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(figsize=(6,6))  # Square figure size

# Display image with value limits set between -0.1 and 0.1
cax = ax.imshow(data_smudge[(int(500000/500)):,:], cmap=colour_style, vmin = -0.04, vmax = 0.05)
ax.set_yticks(np.linspace(0, 1000, 6))
ax.set_yticklabels(np.linspace(1000, 2000, 6).astype(int))
#[:(int(500000/500)),:]

# Add colorbar
plt.colorbar(cax)

ax.set_xlabel("Index along polymer")
ax.set_ylabel(f"Time step")
#ax.set_title(f"Averaged twist along polymer; bracket size: {1+(2*smudge)}")

# Save and show plot
plt.savefig(f"noAxes/kymosmudge{smudge}_{plot_label}_{colour_style}.png", dpi=500, bbox_inches='tight')
plt.show()
