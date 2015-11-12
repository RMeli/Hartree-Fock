import numpy as np
import scipy.special as spec
import scipy.misc as misc

def gaussian_product(aa,bb,Ra,Rb):

    Ra = np.asarray(Ra)
    Rb = np.asarray(Rb)

    R = (aa * Ra + bb * Rb) / (aa + bb)

    c = np.dot(Ra-Rb,Ra-Rb)
    c *= - aa*bb / (aa + bb)
    c = np.exp(c)

    return R,c

def norm(ax,ay,az,aa):
    N = (2*aa/np.pi)**(3./4.)
    N *= (4*aa)**((ax+ay+az)/2)
    N /= np.sqrt( misc.factorial2(2*ax-1) * misc.factorial2(2*ay-1) * misc.factorial2(2*az-1) )

    return N

def Sxyz(a,b,aa,bb,Ra,Rb,R):
    """
    Compute overlap integral between two unnormalized Cartesian gaussian functions along one direction
    """
    S = 0

    for i in range(a+1):
        for j in range(b+1):
            if (i+j) % 2 == 0:
                tmp = misc.comb(a,i,exact=True)
                tmp *= misc.comb(b,j,exact=True)
                tmp *= misc.factorial2(i + j - 1,exact=True)
                tmp /= (2.*(aa+bb))**((i+j)/2.)
                tmp *= (R-Ra)**(a-i)
                tmp *= (R-Rb)**(b-j)

                S += tmp

    return S

def overlap(ax,ay,az,bx,by,bz,aa,bb,Ra,Rb):
    """
    Compute overlap integral between two Cartesian gaussian functions

    Source:
        The Mathematica Journal
        Evaluation of Gaussian Molecular Integrals
        I. Overlap Integrals
        Minhhuy Hô and Julio Manuel Hernández-Pérez
    """

    R,c = gaussian_product(aa,bb,Ra,Rb)

    Na = norm(ax,ay,az,aa)
    Nb = norm(bx,by,bz,bb)

    S = 1
    S *= Sxyz(ax,bx,aa,bb,Ra[0],Rb[0],R[0]) # Sx
    S *= Sxyz(ay,by,aa,bb,Ra[1],Rb[1],R[1]) # Sy
    S *= Sxyz(az,bz,aa,bb,Ra[2],Rb[2],R[2]) # Sz
    S *= Na * Nb * c # Eab and normalization
    S *= (np.pi / (aa+bb))**(3./2.)

    return S

def kinetic(ax,ay,az,bx,by,bz,aa,bb,Ra,Rb):
    """
    Compute kinetic integral between two Cartesian gaussian functions

    Source:
        The Mathematica Journal
        Evaluation of Gaussian Molecular Integrals
        I. Kinetic-Energy Integrals
        Minhhuy Hô and Julio Manuel Hernández-Pérez
    """

    R,c = gaussian_product(aa,bb,Ra,Rb)

    def Kxyz(ac,a1,a2,bc,b1,b2,aa,bb,Ra,Rb,Ra1,Rb1,Ra2,Rb2,Rc,R1,R2):
        kc = 0
        kc += ac * bc * Sxyz(ac-1,bc-1,aa,bb,Ra,Rb,Rc)
        kc += -2 * aa * bc * Sxyz(ac+1,bc-1,aa,bb,Ra,Rb,Rc)
        kc += -2 * ac * bb * Sxyz(ac-1,bc+1,aa,bb,Ra,Rb,Rc)
        kc += 4 * aa * bb * Sxyz(ac+1,bc+1,aa,bb,Ra,Rb,Rc)
        kc *= 0.5

        Kc = 1
        Kc *= c * (np.pi / (aa+bb))**(3./2.) * kc
        Kc *= Sxyz(a1,b1,aa,bb,Ra1,Rb1,R1)
        Kc *= Sxyz(a2,b2,aa,bb,Ra2,Rb2,R2)

        return Kc

    Kx = Kxyz(ax,ay,az,bx,by,bz,aa,bb,Ra[0],Rb[0],Ra[1],Rb[1],Ra[2],Rb[2],R[0],R[1],R[2])
    Ky = Kxyz(ay,az,ax,by,bz,bx,aa,bb,Ra[1],Rb[1],Ra[2],Rb[2],Ra[0],Rb[0],R[1],R[2],R[0])
    Kz = Kxyz(az,ax,ay,bz,bx,by,aa,bb,Ra[2],Rb[2],Ra[0],Rb[0],Ra[1],Rb[1],R[2],R[0],R[1])

    Na = norm(ax,ay,az,aa)
    Nb = norm(bx,by,bz,bb)

    K = (Kx + Ky + Kz) * Na * Nb

    return K
