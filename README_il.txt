#!/bin/bash

#########################################################################
##
## Simulating the effects of a torsion on the polymer, emulating the
## action of a gyrase enyzme
## This is achieved by defining a torsion type interaction, which
## is switched on partway through the simulation
## 
## Initial configuration already provided
##
#########################################################################


#########################################################################
# 1. Run a set of lammps simulations

# The lammps script loops 10 times over, with 10 randomised seeds
# NB: the seeds were input manually so there are only 10; running 11 or 
# more will result in identical results

# LAMMPS is run in 'partition' mode, where 10 simulations run at the
# same time, with a different value of bending stiffness is used in each
# simulation.

# Before running LAMMPS, need to make sure no conda environment running
# and create directories for all the outputs

conda deactivate

mkdir dumpfiles restarts

mpirun -n 7 ~/lammps/src/lmp_mpi -p 7x1 -in in.lammps_SCloop

# Note: Within the lammps script in.lammps_SCloop there is a
#       loop which runs 10 repeat simulations, one after the other
#       using a different seed for the random numbers.
#       Since 11 simulations run at the same time, each with different
#       Lk, then this command in total runs 70 simulations

# This LAMMPS script runs for 30 million steps.
# Scaling based on past simulations, we might expect this to take 
# around 55 minutes per simulation, and 10 hours in total

# You definitely don't want to sit and wait for this to run.
# You can instead use the command in this form so that it runs in the
# background and you can log out and come back later:

nohup mpirun -n 7 ~/newLammps/lmp_mpi -p 7x1 -in in.lammps_SCloop &

# see the page on Teams about 'Some useful command line stuff' for 
# detail on nohup ... &
#########################################################################

#########################################################################
# 2. Process dumpfiles, getting the values for the entire polymer,
     and data for twist in individual sections

# Again, make sure no conda environment is running and create directory:

conda deactivate

mkdir measurements split_measurements

# Loop over torque values and runs

for torque in 1 2 3 4 5 10

do 
    for run in {1..10}
    do
        python process_dump_file.py dumpfiles/dump_L100_t${torque}_run${run}.DNA 100 5000 measurements/scMeasure_L100_t${torque}_run${run}.dat
	python process_dump_file_partition.py dumpfiles/dump_L100_t${torque}_run${run}.DNA 100 5 5000 partition_measurements/scMeasure_L100_t${torque}_run${run}.dat
    done
done

# Note: Again the whole set of commands with two nested for loops
#       was pasted onto the command line. So the process dump code
#       was run 70 times, once for each dump file, one after the other
#########################################################################

#########################################################################
4. Plot single set of runs (8 simulations) to see the long term
   behaviour of the parameters: twist and rg2 for both the stiff and
   normal sections

python plots_vs_time.py

#########################################################################

#########################################################################
# 5. Process the data in full. This version of the processor outputs a
     large selection of data:
	- timestep
	- The rg2 of the normal section
	- The rg2 of the stiff section
	- Total Twist
	- Twist in the normal section
	- Twist in the stiff section
	- The "twist deficit", which is the twist contained in the
	  boundary between the normal and stiff sections, which
	  has to be accounted for to reconcile the previous two
	  series with the total twist
	- Total writhe
	- Writhe in the stiff section
	- Writhe in the normal section
     NB: the writhe in the sections was calculated using an
     approximation in which atoms are near-linearly interpolated
     between the endpoints of a given section.

# Again, make sure no conda environment is running and create directory:

conda deactivate

mkdir split_measurements

# Loop over bending stiffness (bs) values and runs

for porsion in 20.0 30.0 40.0 50.0 60.0 70.0 80.0

do
    for run in {1..10}
    do
        python process_dump_file_split.py dumpfiles/dump_L300_porsion${porsion}_run${run}.DNA 300 150 50000 split_measurements/scMeasure_L300_porsion${porsion}_run${run}.dat
    done
done

# Note: The extra argument specifies how many atoms were included in the 
# bond altering. For this simulation the value was 150.

# This data analysis took a while, around 5 hours to run in total. It 
# might be worth looking at implementing an mpi process to streamline 
# this
#########################################################################

python process_dump_file_split_partition.py dumpfiles/dump_L300_porsion80.0_run5.DNA 300 150 50 50000 partition_measurements/partition_data80.dat

#########################################################################
6. Generate histograms for the twist and rg2 for both the twisted
   and untwisted sections

python plots_histograms_stats.py

# This was modified to reflect the different data we wish to extract 


