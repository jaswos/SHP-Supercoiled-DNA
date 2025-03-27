#!/usr/bin/env python
"""
Program to read in a dump file frame by frame
and do some calculations on what is assumed to be an
inhomogeneous polymer. The writhe data has not been 
as promising as the twist for all series but
was most convincing for the high bending 
stiffness series

Assumes the dump file has the following format
id type xs ys zs ix iy iz c_quat[1] c_quat[2] c_quat[3] c_quat[4]
were xs etc. are in units where the range from 0 to 1
need to adjust these to range from -L/2 to L/2

This code outputs the following values for each timestep (only
a few are actually useful, the rest mostly come in handy when 
checking that the linking number is conserved):
    
{NRg} {SRg} {ATw} {NTw} {STw} {Tw_deficit} {AWr} {SWr} {NWr}

- NRg/ SRg; the normal and stiff radiuses of gyration respectively
- ATw/NTw/STw; the total (all), normal, and stiff twists. Note for the 
twist to be calculated in open sections a modified version of twist.py
was used
- Tw_deficit; the twist "lost" when considering the two separate halves,
specifically, the twist in the joints between the two sections
- AWr/NWr/SWr; the total, normal, and stiff writhes. Note this required the
development of an approximation where a noisy interpolation is added between
the endpoints of the sections; writhe is a property of a closed loop
so passing in an open section wouldn't work.

The twist code was made following the algorith described in "Characterization 
of the Geometry and Topology of DNA Pictured As a Discrete Collection of 
Atoms", Clauvelin (2012).
"""

import sys
import numpy as np
import operator
import os
import random

from supercoiling.twist import calculate_Twist_loop
from supercoiling.twist import calculate_Twist_section
from supercoiling.writhe_cython import calculate_Writhe


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



def readframe_unwrap(infile, N, Naltered):
    """ Read a single frame of atoms from a dump file 
        Rescale coordinates to be in range -L/2 to L/2
        Unwrap coordinates for periodic box """


    """
    This selection of arrays and their construction is very specific to the 
    methods used in the paper Clauvelin et al., J. Chem. TheoryComput. 8 (2012) 1092-1107.
    In this method, three consecutive atoms are needed to calculate a twist
    in the polymer. Hence, the individual sections need an additional atom at the start
    that isn't part of the stiff/normal section, but contributes in calculating the twist .
    Similarly, although only two atoms constitute the boundary between the stiff/normal
    sections, 3 atoms are needed to calculate the vectors necessary.
    """
    Aatoms = [] # all atoms; this will be used to compare total twist/writhe
    Natoms = [] # normal atoms; ones for which id is greater than Naltered
    Satoms = ["placeholder"] # stiffer atoms; ones up to Naltered
    middle = []
    end = [N-1, N, 1]
    # unlike the other arrays, end passe over the end of the loop, so
    # has to be formatted/sorted differently
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
            
            
        #see the docstring above for rationale of these selections 
        if newatom.id == N: 
            Satoms[0] = newatom
        
        if newatom.id == Naltered:
            Natoms.append(newatom)
            
        if newatom.id < Naltered + 1:    
            Satoms.append(newatom)
        else:
            Natoms.append(newatom)
            
        # Need 3 atoms to calculate twist, so take the 3 at the end
        # and at the boundary for this 
        if newatom.id == N-1:
            end[0] = newatom
        if newatom.id == N:
            end[1] = newatom
        if newatom.id == 1:
            end[2] = newatom
        
        if newatom.id == Naltered-1 or newatom.id == Naltered or newatom.id == Naltered + 1:
            middle.append(newatom)
            
            
        Aatoms.append(newatom) 
        # this array should end up identical to how the original data would
        # have parsed the file

    # make sure atoms are sorted by id
    Aatoms.sort(key=operator.attrgetter('id'))
    Natoms.sort(key=operator.attrgetter('id'))
    
    """
    Satoms requires a little rearrangements; atom N is sorted the end but the reason
    for inserting it is to allow the computation of angle 1. We shift atom N to the 0th
    position, and the rest of the array gets shifted up
    """
    Satoms.sort(key=operator.attrgetter('id'))
    last_element = [Satoms.pop()]
    Satoms = last_element + Satoms
    
    middle.sort(key=operator.attrgetter('id'))

    return Aatoms, Natoms, Satoms, middle, end, L


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


def interpolate_atoms(section):
    """

    Parameters
    ----------
    a : array of objects, class atom.
        This code creates additional atoms in between the endpoints of this
        array to enable the writhe code to be run on it
    Returns an array of atoms, where atom 1 and atom N are now "bonded" 
    enough to passed into the writhe function anyway.
    """
    
    start_pos = section[0].x
    end_pos = section[-1].x
    joining_vector = end_pos - start_pos
    unit_join = joining_vector/(np.linalg.norm(joining_vector))
    
    interpolation = []
    N_interpolated = int(np.linalg.norm(joining_vector))
    
    for index in range(1,N_interpolated+1):
        newatom = Atom()
        
        noise = [0, 0, 0]
        while noise[0] == 0 or noise[1] == 0 or noise[2] ==0:
            for i in range(3):    
                noise[i] = random.uniform(-0.03, 0.03)
        
        # the only parameter needed by the writhe function is the position
        newatom.x = (start_pos + index * unit_join) + noise
        newatom.id = 300 + index
        
        interpolation.append(newatom)

    return interpolation

def calculate_writhe_section(section, interpolation, reverse):
    """

    Parameters
    ----------
    section : array
        This is an open section of polymer, which will have an interpolation
        added to it accordingly
    interpolation : array
        These atoms stretch from end to end and are appended to the section
        according to which half is being treated; depending on which section
        the interpolation is claculated from, these atoms can either simply be
        added in order, or in reverse order (to reflect the opposite ordering
        of the polymer).
    reverse : Boolean
        Tells the code whether to reverse the order in which interpolated atoms
        are added to the seciton

    Returns
    -------
    Writhe for the section

    """
    
    if reverse == True:
        for i in range(len(interpolation)):
            section.append(interpolation[-i])
    elif reverse == False:
        for i in range(len(interpolation)):
            section.append(interpolation[i])
    
    writhe = calculate_Writhe(section)
    
    #for i in range(len(section)):
        #print(f"{section[i].id} to 300: {np.linalg.norm(section[i].x - section[149].x)}")
    
    return writhe
        

def calculate_fvector(atoms):
    
    fvector = []
    for i in range(len(atoms)):
        fvector.append(atoms[i].quat2faxis())
        
    return fvector
        
    

############################################################################
############################################################################
############################################################################

#### Start of the main program

#### Read some things from the command line
## (does this in a dumb way with little error checking)
#### assume order is dumpfile Natoms outfileStart

if not len(sys.argv) == 6:
    print("Command usage: ./process_dump_file_split.py dumpfile Natoms Naltered dumpInterval outfile")
    print("           where   dumpfile       is a file name")
    print("                   Natoms         is an integer")
    print("                   Naltered       is an integer")
    print("                   dumpInterval   is an integer, number of steps between dumps")
    print("                   outfile        is a file name which must not exist")
    exit()

dumpfilename = sys.argv[1]
N = int(sys.argv[2])
Naltered = int(sys.argv[3])
dumpInterval = int(sys.argv[4])
outfile = sys.argv[5]

if not os.path.isfile(dumpfilename):
    print("Cannot find file %s.  Exiting..."%dumpfilename)
    exit()
    
#if os.path.isfile(outfile):
    #print("File %s already exists, will not overwrite. Exiting..."%outfile)
    #exit()
    

Nlines = lines_in_file(dumpfilename)  # get length of file
Nframes = int( Nlines / (N+9) )  # there are 9 header lines in each frame


# open the intput file
inf = open(dumpfilename, 'r')  

# open the output files and print a header
ouf = open(outfile, 'w')  
ouf.write("#, timestep, radius_of_gyration_normal, radius_of_gyration_stiff, Twist_total, Twist_normal, Twist_stiff, Twist_deficit, Writhe_total, Writhe_stiff, Writhe_normal\n")


# go through the file frame by frame
for frame in range(Nframes):
    """
    From here, the following letters will stand for the following things:
        A - All (the entire polymer)
        N - Normal, the half of the polymer that is not stiff
        S - Stiff, the half that is
        m - the 2 middle atoms that bridge the normal and stiff sections
        e - end, the two atoms at either end
    """
    # read the frame, unwrapping periodic coordinates
    Aatoms, Natoms, Satoms, middle, end, L = readframe_unwrap(inf,  N, Naltered)

    # calculate radius of gyration
    ARg = radius_of_gyration(Aatoms, L)
    SRg = radius_of_gyration(Satoms,L)
    NRg = radius_of_gyration(Natoms, L)

    Afvector = calculate_fvector(Aatoms)
    Sfvector = calculate_fvector(Satoms)
    Nfvector = calculate_fvector(Natoms)
    mfvector = calculate_fvector(middle)
    efvector = calculate_fvector(end)

    
    # calculate twist and writhe
    ATw = calculate_Twist_loop(Aatoms, Afvector)
    STw = calculate_Twist_section(Satoms, Sfvector)
    NTw = calculate_Twist_section(Natoms, Nfvector)
    
    # this is twist contained in the between the normal and stiff sections
    mTw = calculate_Twist_section(middle, mfvector)
    eTw = calculate_Twist_section(end, efvector)
    Tw_deficit = mTw + eTw


    # the interpolation is calculated from the stiff atoms section, which means
    # it will be ordered from atom index N to atom index N_altered
    interpolation = interpolate_atoms(Satoms)
    
    AWr = calculate_Writhe(Aatoms)
    #print(f"{AWr=}")
    
    SWr = calculate_writhe_section(Satoms, interpolation, True)
    NWr = calculate_writhe_section(Natoms, interpolation, False)
    #print(f"{SWr=}   {NWr=}")

    #output
    ouf.write( f"{frame*dumpInterval} {NRg} {SRg} {ATw} {NTw} {STw} {Tw_deficit} {AWr} {SWr} {NWr}\n")

# close the files
inf.close()

ouf.close()

# Finished!
