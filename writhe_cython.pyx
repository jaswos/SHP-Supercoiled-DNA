import numpy as np
import sys
cimport cython
cimport numpy as np
cimport libc.math
cimport libc.stdio
cimport libc.stdlib

ctypedef np.float64_t DOUBLE

cdef DOUBLE SMALL = 1e-10

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void make_unit(DOUBLE a[3]):
    """ Convert the vector a into a unit vector with 
    same direction """

    cdef:
        DOUBLE L
        unsigned int i

    L = libc.math.sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] )
    for i in range(3):
        a[i] /= L



@cython.boundscheck(False)
@cython.wraparound(False)
cdef DOUBLE fix_roundoff11(DOUBLE a):
    """ checks for round off in a number in the interval [-1,1]
    before using arcsin or arccos """

    if libc.math.fabs(a)>1.001:
        libc.stdio.printf("Error - number should be in interval [-1,1]")
        libc.stdlib.exit(0)
  
    if a>1: 
        a = 1.0
    if a<-1: 
        a = -1.0

    return a


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void cross(DOUBLE a[3],DOUBLE b[3],DOUBLE c[3]):
    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = a[2]*b[0] - a[0]*b[2]
    c[2] = a[0]*b[1] - a[1]*b[0]


@cython.boundscheck(False)
@cython.wraparound(False)
cdef DOUBLE dot(DOUBLE a[3],DOUBLE b[3]):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

@cython.boundscheck(False)
@cython.wraparound(False)
cdef DOUBLE length(DOUBLE a[3]):
    return libc.math.sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] )

@cython.boundscheck(False)
@cython.wraparound(False)
cdef DOUBLE sign(DOUBLE a):
    if a>0:
        return 1.0
    else:
        return -1.0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calculate_Writhe(a):
    """ function to calculate the writhe of a loop from the 
    atom positions 
    This comes from the paper :
    Klenin and Langowski, Biopolymers. 54 (2000) p307 """

    cdef:
        unsigned int N, i, j, k
        DOUBLE Wr, omega, n1n2, n2n3, n3n4, n4n1
        DOUBLE one[3]
        DOUBLE two[3]
        DOUBLE three[3]
        DOUBLE four[3]
        DOUBLE r12[3]
        DOUBLE r23[3]
        DOUBLE r13[3]
        DOUBLE r34[3]
        DOUBLE r24[3]
        DOUBLE r14[3]
        DOUBLE n1[3]
        DOUBLE n2[3]
        DOUBLE n3[3]
        DOUBLE n4[3]
        DOUBLE cvec[3]


    N = len(a)
    Wr = 0.0

    # get an array of the atoms
    atom = np.zeros((N,3),dtype=np.float64)
    for i in range(N):
        atom[i] = a[i].x
    cdef DOUBLE[:,:] at = atom


    for i in range(1,N):
        for  j in range (0,i):

            if (i==j) or (j==i-1) or (i==N-1 and j==0):
                omega = 0.0
            else:
                if j==0:
                    for k in range(3):
                        three[k] = at[N-1,k]
                else:
                    for k in range(3):
                        three[k] = at[j-1,k]

                for k in range(3):
                    one[k] = at[i-1,k]
                    two[k] = at[i,k]
                    four[k] = at[j,k]

                    r12[k] = two[k]-one[k]
                    r34[k] = four[k]-three[k]
                    r23[k] = three[k]-two[k]
                    r13[k] = three[k]-one[k]
                    r14[k] = four[k]-one[k]
                    r24[k] = four[k]-two[k]

                cross(r13,r14,n1) 
                if length(n1)<SMALL:
                    libc.stdio.printf("Error in writhe 1 %.5f",length(n1))
                make_unit(n1)

                cross(r14,r24,n2)
                if length(n2)<SMALL:
                    libc.stdio.printf("Error in writhe 2 %.5f",length(n2))
                make_unit(n2)

                cross(r24,r23,n3)
                if length(n3)<SMALL:
                    libc.stdio.printf("Error in writhe 3 %.5f",length(n3))
                make_unit(n3)

                cross(r23,r13,n4)
                if length(n4)<SMALL:
                    libc.stdio.printf("Error in writhe 4 %.5f",length(n4))
                make_unit(n4)

                n1n2 =  fix_roundoff11( dot(n1,n2) )
                n2n3 =  fix_roundoff11( dot(n2,n3) )
                n3n4 =  fix_roundoff11( dot(n3,n4) )
                n4n1 =  fix_roundoff11( dot(n4,n1) )

                cross(r34,r12,cvec)
                
                omega = ( libc.math.asin( n1n2 ) + libc.math.asin( n2n3 ) + 
                          libc.math.asin( n3n4 ) + 
                          libc.math.asin( n4n1 ) )*sign( dot(cvec,r13) ) 
	

            Wr+= omega/(4.0*libc.math.M_PI)


    Wr*=2.0

    return Wr
