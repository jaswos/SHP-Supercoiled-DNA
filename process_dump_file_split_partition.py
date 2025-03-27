#!/usr/bin/env python

# Program to read in a dump file frame by frame
# and do some calculations.

# Assumes the dump file has the following format
# id type xs ys zs ix iy iz c_quat[1] c_quat[2] c_quat[3] c_quat[4]
# were xs etc. are in units where the range from 0 to 1
# need to adjust these to range from -L/2 to L/2


import sys
import numpy as np
import operator
import os

from supercoiling.twist import calculate_Twist_loop
from supercoiling.twist import calculate_Twist_section
#from supercoiling.writhe_cython import calculate_Writhe


class Atom:
    """ A Class for storing atom information """

    def __init__(self):
        """ Initialise the class """
        self.id = 0                              # id of the atom
        self.type = 0                            # type of the atom
        self.x = np.array([0.0,0.0,0.0],dtype=np.float64)     # position of the atom
        self.image = np.array([0,0,0],dtype=np.int32)         # image flags for atoms
        self.q = np.array([1.0,0.0,0.0,0.0],dtype=np.float64) # quaternion of the atom


    def sep(self,B):
        """ Find separation of this Atom and Atom B """
        return np.sqrt( (self.x[0]-B.x[0])**2 + 
                        (self.x[1]-B.x[1])**2 + (self.x[2]-B.x[2])**2 )

    def minus(self,B):
        """ Subtract B.x from self.x vectorially """
        AminusB = np.array([0.0,0.0,0.0],dtype=np.float64)
        for i in range(3):
            AminusB[i] = self.x[i] - B.x[i]
        return AminusB

    def xdot(self,B):
        """ Find dot product of position x of this Atom and Atom B """
        AdotB = np.array([0.0,0.0,0.0],dtype=np.float64)
        # TO DO : find AdotB
        return AdotB

    def quat2uaxis(self):
        """ Return the unit vector corresponding to the z-axis of the bead """
        p = np.array([0.0,0.0,0.0],dtype=np.float64)
        p[0] = 2.0*self.q[1]*self.q[3] + 2.0*self.q[0]*self.q[2]
        p[1] = 2.0*self.q[2]*self.q[3] - 2.0*self.q[0]*self.q[1]
        p[2] = self.q[0]*self.q[0] - self.q[1]*self.q[1] - self.q[2]*self.q[2] + self.q[3]*self.q[3]
        return p

    def quat2faxis(self):
        """ Return the unit vector corresponding to the z-axis of the bead """
        p = np.array([0.0,0.0,0.0],dtype=np.float64)
        p[0] = self.q[0]*self.q[0] + self.q[1]*self.q[1] - self.q[2]*self.q[2] - self.q[3]*self.q[3]
        p[1] = 2.0*self.q[1]*self.q[2] + 2.0*self.q[0]*self.q[3]
        p[2] = 2.0*self.q[1]*self.q[3] - 2.0*self.q[0]*self.q[2]
        return p

    def quat2vaxis(self):
        """ Return the unit vector corresponding to the z-axis of the bead """
        p = np.array([0.0,0.0,0.0],dtype=np.float64)
        p[0] = 2.0*self.q[1]*self.q[2] - 2.0*self.q[0]*self.q[3]
        p[1] = self.q[0]*self.q[0] - self.q[1]*self.q[1] + self.q[2]*self.q[2] - self.q[3]*self.q[3]
        p[2] = 2.0*self.q[2]*self.q[3] + 2.0*self.q[0]*self.q[1]
        return p



def readframe_unwrap(infile, N):
    """ Read a single frame of atoms from a dump file 
        Rescale coordinates to be in range -L/2 to L/2
        Unwrap coordinates for periodic box """

    atoms = []
    L = []

    # read in the 9 header lines of the dump file
    for i in range(9):
        line = infile.readline()
        if i==5 or i==6 or i==7:
            # get the box size
            line = line.split()
            L.append( np.float64(line[1]) - np.float64(line[0]) )

    # now read the atoms
    for i in range(N):
        line = infile.readline()
        line = line.split()
        newatom = Atom()
        newatom.id = int(line[0])
        newatom.type = int(line[1])
        for j in range(3):
            newatom.x[j] = (np.float64(line[j+2])-0.5)*L[j] # scale
        for j in range(3):
            newatom.image[j] = np.int32(line[j+5])
            newatom.x[j] = newatom.x[j] + newatom.image[j]*L[j] # unwrap
        for j in range(4):
            newatom.q[j] = np.float64(line[j+8])
            
        atoms.append(newatom) 
        # this array should end up identical to how the original data would
        # have parsed the file

    # make sure atoms are sorted by id
    atoms.sort(key=operator.attrgetter('id'))

    return atoms, L


def lines_in_file(filename):
    """ Get the number of lines in the file """

    with open(filename) as f:
        for i, l in enumerate(f):
            pass

    return i + 1


def radius_of_gyration(atoms,L):
    """ Calculate the radius of gytation -- Rg^2 = (1/N) sum ( r_k - r_mean)^2 
    remember to unwrap periodic boundaries """

    # get mean position
    r_mean = np.zeros(3,dtype=np.float64)
    for i in range(len(atoms)):
        r_mean[0] += atoms[i].x[0]
        r_mean[1] += atoms[i].x[1]
        r_mean[2] += atoms[i].x[2]
    r_mean = r_mean/len(atoms)

    # get Rg2
    Rg2 = 0.0
    for i in range(len(atoms)):
        Rg2 += np.sum( np.square( atoms[i].x - r_mean ) )
    Rg2 = Rg2/len(atoms)

    return np.sqrt( Rg2 )

def calculate_partition_twist(atoms, fvector, partition_indices, N_partitions):
    
    # each partition's total twist will be asigned here
    partition_twists = np.zeros(partition_indices.shape[1])

    # partition index array used to create arrays with the necessary
    # atoms for each section, from which twist is calculated
    for i in range(N_partitions):
        atom_selection = []
        fvector_selection = []
        for j in range(partition_indices.shape[0]):
            # somehow the array entries became floats along the journey
            atom_selection.append(atoms[int(partition_indices[j, i])])
            fvector_selection.append(fvector[int(partition_indices[j, i])])
        partition_twists[i] = calculate_Twist_section(atom_selection, fvector_selection)

    return partition_twists    

def calculate_partition_indices(N, partition):
    """

    Parameters
    ----------
    N : total number of atoms 
    partition : the number of atoms in each partition

    Returns
    -------
    array with a number of columns equal to the number of partitions
    where the contents of each column contain all the atom indices 
    needed to calculate the angles for that section

    """
    
    N_partitions = int(N/partition)
    
    # create an array which contains all the indices needed
    # for each partition
    
    indices_needed = np.zeros((partition + 2, N_partitions))
    for i in range(0, N_partitions):
        indices_needed[0, i] = (i * partition) - 1
        
        # this "reverse modulo" acocunts for the fact that index -1 is really
        # referring to the atom at the end of the chain
        if indices_needed[0, i] == -1:
            indices_needed[0, i] = N-1
            
        for j in range(partition):
            indices_needed[j+1, i] = ((i * partition) + j)
        indices_needed[-1, i] = ((i * partition) + partition)
        
        # again, index 300 really refers to index 0; in this
        # sense we're working in a regular mod 300 system
        if indices_needed[-1, i] == N:
            indices_needed[-1, i] = 0
    
    return indices_needed
    

############################################################################
############################################################################
############################################################################

#### Start of the main program

#### Read some things from the command line
## (does this in a dumb way with little error checking)
#### assume order is dumpfile Natoms outfileStart

if not len(sys.argv) == 7:
    print("Command usage: ./process_dump_file.py dumpfile Natoms dumpInterval outfile")
    print("           where   dumpfile       is a file name")
    print("                   Natoms         is an integer")
    print("                   Naltered       is an integer")
    print("                   partition      is an integer")
    print("                   dumpInterval   is an integer, number of steps between dumps")
    print("                   outfile        is a file name which must not exist")
    exit()

dumpfilename = sys.argv[1]
N = int(sys.argv[2])
Naltered = int(sys.argv[3])
partition = int(sys.argv[4])
dumpInterval = int(sys.argv[5])
outfile = sys.argv[6]

if not os.path.isfile(dumpfilename):
    print("Cannot find file %s.  Exiting..."%dumpfilename)
    exit()

if N % partition != 0:
    print(f"Choose partition which is a factor of {N} and, if possible, of {Naltered} as well")
    exit()
    
#if os.path.isfile(outfile):
    #print("File %s already exists, will not overwrite. Exiting..."%outfile)
    #exit()
    

Nlines = lines_in_file(dumpfilename)  # get length of file
Nframes = int( Nlines / (N+9) )  # there are 9 header lines in each frame
N_partitions = int(N/partition) # 

# open the intput file
inf = open(dumpfilename, 'r')  

# open the output files and print a header
ouf = open(outfile, 'w')  

heading = "timestep, Total_Twist"
for part in range(N_partitions):
    heading += f", partition_{part+1}_Tw"
heading += "\n"
ouf.write(heading)



partition_indices = calculate_partition_indices(N, partition)

# go through the file frame by frame
for frame in range(Nframes):
    
    # read the frame, unwrapping periodic coordinates
    atoms, L = readframe_unwrap(inf, N)

    fvector = []
    for i in range(len(atoms)):
        fvector.append(atoms[i].quat2faxis())
    
    # calculate twist in total and in all partitions
    Tw = calculate_Twist_loop(atoms, fvector)
    partition_twist = calculate_partition_twist(atoms, fvector, partition_indices, N_partitions)

    # output
    output = f"{frame*dumpInterval} {Tw}"
    for part in range(N_partitions):
        output += f" {partition_twist[part]}"
    output += "\n"
    ouf.write(output)


# close the files
inf.close()

ouf.close()

# Finished!
