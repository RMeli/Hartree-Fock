"""
    Copyright (C) 2015 Rocco Meli

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import scipy.misc as misc
import scipy.special as spec
import scipy.integrate as quad

from basis import *

def gaussian_product(aa,bb,Ra,Rb):
    """
    Gaussian produc theorem.

    INPUT:
        AA: Exponential coefficient of Gaussian 1
        BB: Exponential coefficient of Gaussian 2
        RA: Center of Gaussian 1
        RB: Center of Gaussian 2
    OUTPUT:
        R: Gaussian product center
        C: Gaussian product coefficient

    Source:
        Modern Quantum Chemistry
        Szabo and Ostlund
        Dover
        1989
    """

    # Transform centers in Numpy arrays
    Ra = np.asarray(Ra)
    Rb = np.asarray(Rb)

    # Compute Gaussian product center
    R = (aa * Ra + bb * Rb) / (aa + bb)

    # Compute Gaussian product coefficient
    c = np.dot(Ra-Rb,Ra-Rb)
    c *= - aa*bb / (aa + bb)
    c = np.exp(c)

    return R,c

def norm(ax,ay,az,aa):
    """
    General cartesian Gaussian normalization factor.

    INPUT:
        AX: Angular momentum lx
        AY: Angular momentum ly
        AZ: Angular momentum lz
        AA: Gaussian exponential coefficient
    OUTPUT:
        N: Normalization coefficient (to be multiplied with the Gaussian)

    Source:
        Handbook of Computational Chemistry
        David Cook
        Oxford University Press
        1998
    """

    # Compute normalization coefficient
    N = (2*aa/np.pi)**(3./4.)
    N *= (4*aa)**((ax+ay+az)/2)
    N /= np.sqrt( misc.factorial2(2*ax-1) * misc.factorial2(2*ay-1) * misc.factorial2(2*az-1) )

    return N

def Sxyz(a,b,aa,bb,Ra,Rb,R):
    """
    Compute overlap integral between two unnormalized Cartesian gaussian functions along one direction.

    INPUT:
        A: Angular momentum along the chosen direction for Gaussian 1
        B: Angular momentum along the chosen direction for Gaussian 2
        AA: Exponential coefficient for Gaussian 1
        BB: Exponential coefficient for Gaussian 2
        RA: Coordinate (along chosen direction) of the center of Gaussian 1
        RB: Coordinate (along chosen direction) of the center of Gaussian 2
        R: Coordinate (along chosen direction) of the center of the product of the two gaussians
    OUTPUT:
        S: Overlap of the two gaussians along the chosen direction

    Source:
        The Mathematica Journal
        Evaluation of Gaussian Molecular Integrals
        I. Overlap Integrals
        Minhhuy Hô and Julio Manuel Hernández-Pérez
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
    Compute overlap integral between two Cartesian gaussian functions.

    INPUT:
        AX: Angular momentum lx for Gaussian 1
        AY: Angular momentum ly for Gaussian 1
        AZ: Angular momentum lz for Gaussian 1
        AA: Exponential coefficient for Gaussian 1
        BX: Angular momentum lx for Gaussian 2
        BY: Angular momentum ly for Gaussian 2
        BZ: Angular momentum lz for Gaussian 2
        BB: Exponential coefficient for Gaussian 2
        RA: Center of Gaussian 1
        RB: Center of Gaussian 2
    OUTPUT:
        S: Overlap of the two gaussians

    Source:
        The Mathematica Journal
        Evaluation of Gaussian Molecular Integrals
        I. Overlap Integrals
        Minhhuy Hô and Julio Manuel Hernández-Pérez
    """

    # Compute gaussian product center and coefficient
    R,c = gaussian_product(aa,bb,Ra,Rb)

    # Compute normalization factors for the two gaussians
    Na = norm(ax,ay,az,aa)
    Nb = norm(bx,by,bz,bb)

    S = 1
    S *= Sxyz(ax,bx,aa,bb,Ra[0],Rb[0],R[0]) # Overlap along x
    S *= Sxyz(ay,by,aa,bb,Ra[1],Rb[1],R[1]) # Overlap along y
    S *= Sxyz(az,bz,aa,bb,Ra[2],Rb[2],R[2]) # Overlap along z
    S *= Na * Nb * c # Product coefficient and normalization
    S *= (np.pi / (aa+bb))**(3./2.) # Normalization

    return S

def kinetic(ax,ay,az,bx,by,bz,aa,bb,Ra,Rb):
    """
    Compute kinetic integral between two Cartesian gaussian functions.

    INPUT:
        AX: Angular momentum lx for Gaussian 1
        AY: Angular momentum ly for Gaussian 1
        AZ: Angular momentum lz for Gaussian 1
        AA: Exponential coefficient for Gaussian 1
        BX: Angular momentum lx for Gaussian 2
        BY: Angular momentum ly for Gaussian 2
        BZ: Angular momentum lz for Gaussian 2
        BB: Exponential coefficient for Gaussian 2
        RA: Center of Gaussian 1
        RB: Center of Gaussian 2
    OUTPUT:
        K: Kinetic integral between the two gaussians

    Source:
        The Mathematica Journal
        Evaluation of Gaussian Molecular Integrals
        II. Kinetic-Energy Integrals
        Minhhuy Hô and Julio Manuel Hernández-Pérez
    """

    R,c = gaussian_product(aa,bb,Ra,Rb)

    def Kxyz(ac,a1,a2,bc,b1,b2,aa,bb,Ra,Rb,Ra1,Rb1,Ra2,Rb2,Rc,R1,R2):
        """
            Compute kinetic integral between two Cartesian gaussian functions along one direction.

            INPUT:
                AC: Component of angular momentum for Gaussian 1 along direction of interest
                A1: Component of angular momentum for Gaussian 1 along second direction
                A2: Component of angular momentum for Gaussian 1 along third direction
                BC: Component of angular momentum for Gaussian 2 along direction of interest
                B1: Component of angular momentum for Gaussian 2 along second direction
                B2: Component of angular momentum for Gaussian 2 along third direction
                AA: Exponential coefficient for Gaussian 1
                BB: Exponential coefficient for Gaussian 2
                RA: Component of the center of Gaussian 1 along direction of interest
                RB: Component of the center of Gaussian 2 along direction of interest
                RA1: Component of the center of Gaussian 1 along second direction
                RA2: Component of the center of Gaussian 1 along third direction
                RB1: Component of the center of Gaussian 2 along second direction
                RB2: Component of the center of Gaussian 2 along third direction
            OUTPUT:
                KC: Kinetic integral between two gaussians along chosen direction

            Source:
                The Mathematica Journal
                Evaluation of Gaussian Molecular Integrals
                II. Kinetic-Energy Integrals
                Minhhuy Hô and Julio Manuel Hernández-Pérez
        """

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

    # Cyclic permutation of the entries
    Kx = Kxyz(ax,ay,az,bx,by,bz,aa,bb,Ra[0],Rb[0],Ra[1],Rb[1],Ra[2],Rb[2],R[0],R[1],R[2]) # Kinetic integral along x
    Ky = Kxyz(ay,az,ax,by,bz,bx,aa,bb,Ra[1],Rb[1],Ra[2],Rb[2],Ra[0],Rb[0],R[1],R[2],R[0]) # Kinetic integral along y
    Kz = Kxyz(az,ax,ay,bz,bx,by,aa,bb,Ra[2],Rb[2],Ra[0],Rb[0],Ra[1],Rb[1],R[2],R[0],R[1]) # Kinetic integral along z

    Na = norm(ax,ay,az,aa) # Normalization factor for Gaussian 1
    Nb = norm(bx,by,bz,bb) # Normalization factor for Gaussian 2

    K = (Kx + Ky + Kz) * Na * Nb # Normalization of total kinetic energy integral

    return K

def f(j,l,m,a,b):
    """
    Expansion coefficient f.

    Source:
        Handbook of Computational Chemistry
        David Cook
        Oxford University Press
        1998
    """

    f = 0

    for k in range(max(0,j-m),min(j,l)+1):
        tmp = 1
        tmp *= spec.binom(l,k)
        tmp *= spec.binom(m,j-k)
        tmp *= a**(l-k)
        tmp *= b**(m+k-j)

        f += tmp

    return f

def F(nu,x):
    """
    Boys function.

    INPUT:
        NU: Boys function index
        X: Boys function variable

    OUTPUT:
        FF: Value of the Boys function for index NU evaluated at X

    Source:
        Evaluation of the Boys Function using Analytical Relations
        I. I. Guseinov and B. A. Mamedov
        Journal of Mathematical Chemistry
        2006
    """

    if x < 1e-8:
        # Taylor expansion for argument close or equal to zero (avoid division by zero)
        ff =  1 / (2 * nu + 1) - x / (2 * nu + 3)
    else:
        # Evaluate Boys function with incomplete and complete Gamma functions
        ff = 0.5 / x**(nu+0.5) * spec.gamma(nu+0.5)*spec.gammainc(nu+0.5,x)

    return ff

def nuclear(ax,ay,az,bx,by,bz,aa,bb,Ra,Rb,Rn,Zn):
    """
    Compute nuclear-electron interaction integrals.

    INPUT:
    AX,AY,AZ: Angular momentum components for the first Gaussian.
    BX,BY,BZ: Angular momentum components for the second Gaussian.
    AA: Exponential coefficient for the first Gaussian.
    BB: Exponential coefficient for the second Gaussian.
    RA: Center of the first Gaussian.
    RB: Center of the second Gaussian.
    RN: Nuclear coordinates.
    ZN: Nuclear charge.

    Source:

        Handbook of Computational Chemistry
        David Cook
        Oxford University Press
        1998
    """

    Vn = 0

    # Intermediate variable
    g = aa + bb
    eps = 1. / (4 * g)

    Rp,c = gaussian_product(aa,bb,Ra,Rb) # Gaussian product

    def A(l,r,i,l1,l2,Ra,Rb,Rc,Rp):
        """
        Expansion coefficient A.

        Source:
            Handbook of Computational Chemistry
            David Cook
            Oxford University Press
            1998
        """

        A = 1
        A *= (-1)**(l)
        A *= f(l,l1,l2,Rp-Ra,Rp-Rb)
        A *= (-1)**i
        A *= misc.factorial(l,exact=True)
        A *= (Rp-Rc)**(l-2*r-2*i)
        A *= eps**(r+i)
        A /= misc.factorial(r,exact=True)
        A /= misc.factorial(i,exact=True)
        A /= misc.factorial(l-2*r-2*i,exact=True)

        return A

    for l in range(0,ax+bx+1):
        for r in range(0,int(l/2)+1):
            for i in range(0,int((l-2*r)/2)+1):
                Ax = A(l,r,i,ax,bx,Ra[0],Rb[0],Rn[0],Rp[0])

                for m in range(0,ay+by+1):
                    for s in range(0,int(m/2)+1):
                        for j in range(0,int((m-2*s)/2)+1):
                            Ay = A(m,s,j,ay,by,Ra[1],Rb[1],Rn[1],Rp[1])

                            for n in range(0,az+bz+1):
                                for t in range(0,int(n/2)+1):
                                    for k in range(0,int((n-2*t)/2)+1):
                                        Az =  A(n,t,k,az,bz,Ra[2],Rb[2],Rn[2],Rp[2])

                                        nu = l + m + n - 2 * (r + s + t) - (i + j + k) # Index of Boys function

                                        ff = F(nu,g*np.dot(Rp-Rn,Rp-Rn)) # Boys function

                                        Vn += Ax * Ay * Az * ff

    # Compute normalization
    Na = norm(ax,ay,az,aa)
    Nb = norm(bx,by,bz,bb)

    Vn *= - Zn * Na * Nb * c * 2 * np.pi / g

    return Vn


def electronic(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,aa,bb,cc,dd,Ra,Rb,Rc,Rd):
    """
    Compute electron-electron interaction integrals.

    INPUT:
        AX,AY,AZ: Angular momentum components for the first Gaussian.
        BX,BY,BZ: Angular momentum components for the second Gaussian.
        CX,CY,CZ: Angular momentum components for the third Gaussian.
        DX,DY,DZ: Angular momentum components for the fourth Gaussian.
        AA: Exponential coefficient for the first Gaussian.
        BB: Exponential coefficient for the second Gaussian.
        CC: Exponential coefficient for the third Gaussian.
        DD: Exponential coefficient for the fourth Gaussian.
        RA: Center of the first Gaussian.
        RB: Center of the second Gaussian.
        RC: Center of the third Gaussian.
        RD: Center of the fourth Gaussian.
    OUTPUT:
        G: Electron-electron integral

    Source:

        Handbook of Computational Chemistry
        David Cook
        Oxford University Press
        1998

        ERRATA (the original formula is WRONG)!
            http://spider.shef.ac.uk/
    """


    G = 0

    # Intermediate variable
    g1 = aa + bb
    g2 = cc + dd

    # Compute gaussian products
    Rp, c1 = gaussian_product(aa,bb,Ra,Rb)
    Rq, c2 = gaussian_product(cc,dd,Rc,Rd)

    delta = 1 / (4 * g1) + 1 / (4 * g2)

    def theta(l,l1,l2,a,b,r,g):
        """
        Expansion coefficient theta.

        Source:
            Handbook of Computational Chemistry
            David Cook
            Oxford University Press
            1998
        """

        t = 1
        t *= f(l,l1,l2,a,b)
        t *= misc.factorial(l,exact=True)
        t *= g**(r-l)
        t /= misc.factorial(r,exact=True) * misc.factorial(l-2*r,exact=True)

        return t


    def B(l,ll,r,rr,i,l1,l2,Ra,Rb,Rp,g1,l3,l4,Rc,Rd,Rq,g2):
        """
        Expansion coefficient B.

        Source:
            Handbook of Computational Chemistry
            David Cook
            Oxford University Press
            1998
        """

        b = 1
        b *= (-1)**(l) * theta(l,l1,l2,Rp-Ra,Rp-Rb,r,g1)
        b *= theta(ll,l3,l4,Rq-Rc,Rq-Rd,rr,g2)
        b *= (-1)**i * (2*delta)**(2*(r+rr))
        b *= misc.factorial(l + ll - 2*r - 2*rr,exact=True)
        b *= delta**i * (Rp-Rq)**(l+ll-2*(r+rr+i))

        tmp = 1
        tmp *= (4*delta)**(l+ll) * misc.factorial(i,exact=True)
        tmp *= misc.factorial(l+ll-2*(r+rr+i),exact=True)

        b /= tmp

        return b

    for l in range(0,ax+bx+1):
        for r in range(0,int(l/2)+1):
            for ll in range(0,cx+dx+1):
                for rr in range(0,int(ll/2)+1):
                    for i in range(0,int((l+ll-2*r-2*rr)/2)+1):
                        Bx = B(l,ll,r,rr,i,ax,bx,Ra[0],Rb[0],Rp[0],g1,cx,dx,Rc[0],Rd[0],Rq[0],g2)

                        for m in range(0,ay+by+1):
                            for s in range(0,int(m/2)+1):
                                for mm in range(0,cy+dy+1):
                                    for ss in range(0,int(mm/2)+1):
                                        for j in range(0,int((m+mm-2*s-2*ss)/2)+1):
                                            By = B(m,mm,s,ss,j,ay,by,Ra[1],Rb[1],Rp[1],g1,cy,dy,Rc[1],Rd[1],Rq[1],g2)

                                            for n in range(0,az+bz+1):
                                                for t in range(0,int(n/2)+1):
                                                    for nn in range(0,cz+dz+1):
                                                        for tt in range(0,int(nn/2)+1):
                                                            for k in range(0,int((n+nn-2*t-2*tt)/2)+1):
                                                                Bz = B(n,nn,t,tt,k,az,bz,Ra[2],Rb[2],Rp[2],g1,cz,dz,Rc[2],Rd[2],Rq[2],g2)

                                                                nu = l+ll+m+mm+n+nn-2*(r+rr+s+ss+t+tt) - (i + j + k)

                                                                ff = F(nu,np.dot(Rp-Rq,Rp-Rq)/(4.*delta))

                                                                G += Bx * By * Bz * ff


    # Compute normalization
    Na = norm(ax,ay,az,aa)
    Nb = norm(bx,by,bz,bb)
    Nc = norm(cx,cy,cz,cc)
    Nd = norm(dx,dy,dz,dd)

    G *= Na * Nb * Nc * Nd * c1 * c2 * 2 * np.pi**2 / (g1 * g2) * np.sqrt(np.pi / (g1 + g2))

    return G


def EE_list(basis):
    """
    Multidimensional array of two-electron integrals.

    INPUT:
        BASIS: Basis set
    OUTPUT:
        EE: list of two-electron integrals, with indices (i,j,k,l)
    """

    # Size of the basis set
    K = basis.K

    # List of basis functions
    B = basis.basis()

    EE = np.zeros((K,K,K,K))

    Nee= 0

    for i,b1 in enumerate(B):
        for j,b2 in enumerate(B):
            for k,b3 in enumerate(B):
                for l,b4 in enumerate(B):

                    Nee += 1

                    # Print update of calculation (can be long)
                    if Nee % 500 == 0:
                        print("     Computed ", Nee, " two-electron integrals of ", K**4, ".",sep='')


                    for a1,d1 in zip(b1["a"],b1["d"]):
                        for a2,d2 in zip(b2["a"],b2["d"]):
                            for a3,d3 in zip(b3["a"],b3["d"]):
                                for a4,d4 in zip(b4["a"],b4["d"]):
                                    # Basis functions centers
                                    R1 = b1["R"]
                                    R2 = b2["R"]
                                    R3 = b3["R"]
                                    R4 = b4["R"]

                                    # Basis functions angular momenta
                                    ax = b1["lx"]
                                    ay = b1["ly"]
                                    az = b1["lz"]

                                    # Basis functions angular momenta
                                    bx = b2["lx"]
                                    by = b2["ly"]
                                    bz = b2["lz"]

                                    # Basis functions angular momenta
                                    cx = b3["lx"]
                                    cy = b3["ly"]
                                    cz = b3["lz"]

                                    # Basis functions angular momenta
                                    dx = b4["lx"]
                                    dy = b4["ly"]
                                    dz = b4["lz"]

                                    tmp = 1
                                    tmp *= d1.conjugate()*d2.conjugate()
                                    tmp *= d3 * d4
                                    tmp *= electronic(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,a1,a2,a3,a4,R1,R2,R3,R4)

                                    EE[i,j,k,l] += tmp

    return EE

def print_EE_list(ee):
    """
    Print list of electron-electron integrals.

    INPUT:
        EE: list of electron-electron integrals (computed by EE_LIST function)
    """

    K = ee.shape[0]

    for i in range(K):
        for j in range(K):
            for k in range(K):
                for l in range(K):
                    print("({0},{1},{2},{3})  {4}".format(i+1,j+1,k+1,l+1,ee[i,j,k,l]))

if __name__ == "__main__":

    """
    Results compared with

        Modern Quantum Chemistry
        Szabo and Ostlund
        Dover
        1989
    """

    # System: HeH+
    HeH = [Atom("He",(0,0,1.4632),2,["1s"]),Atom("H",(0,0,0),1,["1s"])]

    sto3g_HeH = STO3G(HeH) # Create basis set

    ee_HeH = EE_list(sto3g_HeH) # Compute electron-electron integrals for HeH+

    print("######################")
    print("Two electron integrals")
    print("######################")

    print("\n HeH")
    print_EE_list(ee_HeH) # Print electron-electron integrals
