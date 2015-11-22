import numpy as np
import scipy.misc as misc
import scipy.special as spec
import scipy.integrate as quad

from basis import *

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
        II. Kinetic-Energy Integrals
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

def f(j,l,m,a,b):
    """
    Expansion coefficient.

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

    def f(t):
        return t**(2*nu) * np.exp(-x*t**2)

    F = quad.quad(f,0,1,epsabs=1e-12,epsrel=1e-12,limit=500)[0]
    """

    tol = 1e-15

    if nu < 10:
        resid = 1

        ex = np.exp(-x) * 0.5

        s = spec.gamma(nu + 0.5) / spec.gamma(nu + 1.5)

        b = ex * s

        i = 1
        while resid > tol:
            s += spec.gamma(nu + 0.5) * x**i / spec.gamma(nu + i + 1.5)

            bnew = ex * s

            resid = abs(b - bnew)

            b = bnew

            i += 1

        F = b
        """
    return F


def nuclear(ax,ay,az,bx,by,bz,aa,bb,Ra,Rb,Rn,Zn):
    """
    Compute nuclear-electron interaction integrals.


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

    g = aa + bb
    eps = 1. / (4 * g)

    Rp,c = gaussian_product(aa,bb,Ra,Rb)

    def A(l,r,i,l1,l2,Ra,Rb,Rc,Rp):
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

                                        nu = l + m + n - 2 * (r + s + t) - (i + j + k)

                                        ff = F(nu,g*np.dot(Rp-Rn,Rp-Rn))

                                        Vn += Ax * Ay * Az * ff

    Na = norm(ax,ay,az,aa)
    Nb = norm(bx,by,bz,bb)

    Vn *= - Zn * Na * Nb * c * 2 * np.pi / g

    return Vn


def electronic(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,aa,bb,cc,dd,Ra,Rb,Rc,Rd):
    """
    Compute electron-electron interaction integrals.


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

    Source:

        Handbook of Computational Chemistry
        David Cook
        Oxford University Press
        1998
    """


    G = 0

    g1 = aa + bb
    g2 = cc + dd

    Rp, c1 = gaussian_product(aa,bb,Ra,Rb)
    Rq, c2 = gaussian_product(cc,dd,Rc,Rd)

    delta = 1 / (4 * g1) + 1 / (4 * g2)

    def theta(l,l1,l2,a,b,r,g):
        t = 1
        t *= f(l,l1,l2,a,b)
        t *= misc.factorial(l,exact=True)
        t *= g**(r-l)
        t /= misc.factorial(r,exact=True) * misc.factorial(l-2*r,exact=True)

        return t


    def B(l,ll,r,rr,i,l1,l2,Ra,Rb,Rp,g1,l3,l4,Rc,Rd,Rq,g2):
        b = 1
        b *= (-1)**(ll) * theta(l,l1,l2,Rp-Ra,Rp-Rb,r,g1)
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
            for i in range(0,int((l-2*r)/2)+1):
                for ll in range(0,cx+dx+1):
                    for rr in range(0,int(ll/2)+1):
                        Bx = B(l,ll,r,rr,i,ax,bx,Ra[0],Rb[0],Rp[0],g1,cx,dx,Rc[0],Rd[0],Rq[0],g2)

                        for m in range(0,ay+by+1):
                            for s in range(0,int(m/2)+1):
                                for j in range(0,int((m-2*s)/2)+1):
                                    for mm in range(0,cy+dy+1):
                                        for ss in range(0,int(mm/2)+1):
                                            By = B(m,mm,s,ss,j,ay,by,Ra[1],Rb[1],Rp[1],g1,cy,dy,Rc[1],Rd[1],Rq[1],g2)

                                            for n in range(0,az+bz+1):
                                                for t in range(0,int(n/2)+1):
                                                    for k in range(0,int((n-2*t)/2)+1):
                                                        for nn in range(0,cz+dz+1):
                                                            for tt in range(0,int(nn/2)+1):
                                                                Bz = B(n,nn,t,tt,k,az,bz,Ra[2],Rb[2],Rp[2],g1,cz,dz,Rc[2],Rd[2],Rq[2],g2)

                                                                nu = l+ll+m+mm+n+nn-2*(r+rr+s+ss+t+tt) - (i + j + k)

                                                                ff = F(nu,np.dot(Rp-Rq,Rp-Rq)*g1*g2/(g1+g2))

                                                                G += Bx * By * Bz * ff



    Na = norm(ax,ay,az,aa)
    Nb = norm(bx,by,bz,bb)
    Nc = norm(cx,cy,cz,cc)
    Nd = norm(dx,dy,dz,dd)

    G *= Na * Nb * Nc * Nd * c1 * c2 * 2 * np.pi**2 / (g1 * g2) * np.sqrt(np.pi / (g1 + g2))

    return G

"""

def electronic(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,aa,bb,cc,dd,Ra,Rb,Rc,Rd):

    G = 0

    gp = aa + bb
    gq = cc + dd

    Rp, Cp = gaussian_product(aa,bb,Ra,Rb)
    Rq, Cq = gaussian_product(cc,dd,Rc,Rd)

    etap = aa * bb / gp
    etaq = cc * dd / gq

    eta = etap * etaq / (etap + etaq)

    def J(l1,l2,l3,l4,Ra,Rb,Rc,Rd,Rp,Rq):

        def T(i1,i2,o1,o2,r1,aa,bb,Ra,Rb):
            T = 1
            T *= (-1)**(o2+r1)
            T *= misc.factorial(o1+o2,exact=True)
            T /= 4**(i1+i2+r1)
            T /= misc.factorial(i1,exact=True)
            T /= misc.factorial(i2,exact=True)
            T /= misc.factorial(o1,exact=True)
            T /= misc.factorial(o2,exact=True)
            T /= misc.factorial(r1,exact=True)
            T *= aa**(o2-i1-r1)
            T *= bb**(o1-i2-r1)
            T *= etap**(2*(i1+i2)+r1)
            T *= (Ra-Rb)**(o1+o2-2*r1)
            T /= misc.factorial(l1-2*i1-o1,exact=True)
            T /= misc.factorial(l2-2*i2-o2,exact=True)
            T /= misc.factorial(o1+o2-2*r1,exact=True)
            T *= misc.factorial(l3,exact=True)
            T *= misc.factorial(l4,exact=True)
            T /= etaq**(l3+l4)

        J = 1
        J *= (-1)**(l1+l2)
        J *= misc.factorial(l1,exact=True)
        J *= misc.factorial(l2,exact=True)
        J /= etap**(l1+l2)
        J *= misc.factorial(l3,exact=True)
        J *= misc.factorial(l4,exact=True)
        J /= etaq**(l3+l4)

        SUM = 0

        for i1 in range(0,int(l1/2)+1):
            for i2 in range(0,int(l2/2)+1):
                for o1 in range(0,l1-2*i1+1):
                    for o2 in range(0,l2-2*i2+1):
                        for r1 in range(0,int((o1+o2)/2)+1)):
                            T1 = T(i1,i2,o1,o2,r1,aa,bb,Ra,Rb)

                            for i3 in range(0,int(l3/2)+1):
                                for i4 in range(0,int(l4/2)+1):
                                    for o3 in range(0,l3-2*i3+1):
                                        for o4 in range(0,l4-2*i4+1):
                                            for r2 in range(0,int((o3+o4)/2)+1):
                                                T2 = T(i3,i4,o3,o4,r2,cc,dd,Rc,Rd)

                                                mu = l1 + l2 + l3 + l4 - 2*(i1+i2+i3+i4) -(o1+o2+o3+o4)

                                                    for u in range(0,int(mu/2)+1):
                                                        T3 = 1
                                                        T3 *= (-1)**u
                                                        T3 *= misc.factorial(mu,exact=True)
                                                        T3 *= eta**(mu-u)
                                                        T3 *= (Rp-Rq)**(mu-2*u)
                                                        T3 /= 4**u * misc.factorial(u,exact=True)
                                                        T3 /= misc.factorial(mu-2*u,exact=True)

                                                        SUM += T1*T2*T3

        J *= SUM

        return J

    Na = norm(ax,ay,az,aa)
    Nb = norm(bx,by,bz,bb)
    Nc = norm(cx,cy,cz,cc)
    Nd = norm(dx,dy,dz,dd)

    G = 1
    G *= J(ax,bx,cx,dx,Ra[0],Rb[0],Rc[0],Rd[0],Rp[0],Rq[0])
    G *= J(ay,by,cy,dy,Ra[1],Rb[1],Rc[1],Rd[1],Rp[1],Rq[1])
    G *= J(az,bz,cz,dz,Ra[2],Rb[2],Rc[2],Rd[2],Rp[2],Rq[2])
    G *= 2* F
    G *= Na*Nb*Nc*Nd * np.pi**(5./2.) / (gp * gq * np.sqrt(gp+gq)) * Cp * Cq
"""



def EE_list(basis):
    """
    Multidimensional array of two-electron integrals.
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

                    if Nee % 500 == 0:
                        print("     Computed ", Nee, " two-electron integrals of ", K**4, ".",sep='')


                    for a1,d1 in zip(b1["a"],b1["d"]):
                        for a2,d2 in zip(b2["a"],b2["d"]):
                            for a3,d3 in zip(b3["a"],b3["d"]):
                                for a4,d4 in zip(b4["a"],b4["d"]):
                                    R1 = b1["R"]
                                    R2 = b2["R"]
                                    R3 = b3["R"]
                                    R4 = b4["R"]

                                    ax = b1["lx"]
                                    ay = b1["ly"]
                                    az = b1["lz"]

                                    bx = b2["lx"]
                                    by = b2["ly"]
                                    bz = b2["lz"]

                                    cx = b3["lx"]
                                    cy = b3["ly"]
                                    cz = b3["lz"]

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

    HeH = [Atom("He",(0,0,1.4632),2,["1s"]),Atom("H",(0,0,0),1,["1s"])]

    sto3g_HeH = STO3G(HeH)

    ee_HeH = EE_list(sto3g_HeH)

    print("######################")
    print("Two electron integrals")
    print("######################")

    print("\n HeH")
    print_EE_list(ee_HeH)
