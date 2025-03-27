
import sys
import numpy as np

SMALL = 1e-10

def make_unit(a):
    """ Convert the vector a into a unit vector with 
    same direction """

    L = np.linalg.norm(a)
    if L != 0:
        a /= L
    return a



def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def calculate_Twist_loop(a,f):
    #""" function to calculate the twist round a closed loop from 
    #the f vectors. 
    #This comes from the paper : 
    #Clauvelin et al., J. Chem. TheoryComput. 8 (2012) 1092-1107

    # This code has been adapted to reflect that this is an open
    # section of polymer; the endpoints are not linked

    N = len(f)
    Tw = 0

    # get a set of tangent vectors
    t = []
    for i in range(0,N-1):
      t.append( make_unit( a[i+1].x-a[i].x ) )
    t.append( a[0].x-a[N-1].x )
    t[-1] = make_unit(t[N-1])

    # a set of m vectors which are perpendicular to t and f
    m = []
    for i in range(N):
        m.append( make_unit( np.cross(t[i],f[i]) )  )

    # a set of binomial vectors for every joint
    b = []
    b.append( make_unit( np.cross(t[N-1],t[0]) ) )
    for i in range(1,N):
        b.append( make_unit( np.cross(t[i-1],t[i]) ) )

    # and a set of turning angles between every joint
    angA = []
    angA.append( np.arccos( t[N-1].dot(t[0]) ) )
    for i in range(1,N):
        angA.append( np.arccos( t[i-1].dot(t[i]) ) )

    # now find a shifted m vector
    mshift = []
    mshift.append( np.dot(rotation_matrix(b[0],angA[0]),m[N-1]) )
    for i in range(1,N):
        mshift.append( np.dot(rotation_matrix(b[i],angA[i]),m[i-1]) )

    for i in range(N):
        #costheta = np.dot(mshift[i],m[i])
        sintheta = np.dot(t[i], np.cross(mshift[i],m[i]) )
        Tw += np.arcsin(sintheta)

    return Tw/(2.0*np.pi)

def calculate_Twist_section(a,f):
    #""" function to calculate the twist round a loop from 
    #the f,u and v vectors. 
    #This comes from the paper : 
    #Clauvelin et al., J. Chem. TheoryComput. 8 (2012) 1092-1107 """

    N = len(f)
    Tw = 0

    # Get a set of tangent vectors

    t = []
    for i in range(0,N-1):
        unit = make_unit( a[i+1].x-a[i].x )
        t.append(unit)

    # This "open" code will never consider the behaviour of the
    # last atom; it would usually be assumed to be linked to the 
    # first but this is NOT the case
    
    # a set of m vectors which are perpendicular to t and f
    m = []
    for i in range(0, N-1):
        m.append( make_unit( np.cross(t[i],f[i]) )  )

    # a set of binormal vectors for every joint
    b = ["dud"]
    # to mantain the similarity to Clauvelin, the ith component
    # should correspond to the correct vectors, which requires
    # this shift by 1
    
    for i in range(1, N-1):
        b.append( make_unit( np.cross(t[i-1],t[i]) ) )

    # and a set of turning angles between every joint
    angA = ["dud"]
    for i in range(1, N-1):
        angA.append( np.arccos( t[i-1].dot(t[i]) ) )

    # now find a shifted m vector
    
    mshift = []
    for i in range(1, N-1):
        mshift.append( np.dot(rotation_matrix(b[i],angA[i]),m[i-1]) )

    #print(f"len mshift: {len(mshift)}")
    #print(f"len m: {len(m)}")
    #print(f"len t: {len(t)}")

    for i in range(1, N-1):
        #costheta = np.dot(mshift[i],m[i])
        sintheta = np.dot(t[i], np.cross(mshift[i-1],m[i]) )
        Tw += np.arcsin(sintheta)

    return Tw/(2.0*np.pi)
