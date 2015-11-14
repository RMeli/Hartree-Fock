import numpy as np
import scipy.misc as misc
import scipy.integrate as quad

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

    """

    f = 0

    for k in range(max(0,j-m),min(j,l)+1):
        tmp = 1
        tmp *= misc.comb(l,k,exact=True)
        tmp *= misc.comb(m,j-k,exact=True)
        tmp *= a**(l-k)
        tmp *= b**(m+k-j)

        f += tmp

    return f

def F(nu,x,dt=1e-3):

    def f(t):
        return t**(2*nu) * np.exp(-x*t**2)

    F = quad.quad(f,0,1)[0]

    return F

def nuclear(ax,ay,az,bx,by,bz,aa,bb,Ra,Rb,Rn,Zn):
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
                A1 = A(l,r,i,ax,bx,Ra[0],Rb[0],Rn[0],Rp[0])

                for m in range(0,ay+by+1):
                    for s in range(0,int(m/2)+1):
                        for j in range(0,int((m-2*s)/2)+1):
                            A2 = A(m,s,j,ay,by,Ra[1],Rb[1],Rn[1],Rp[1])

                            for n in range(0,az+bz+1):
                                for t in range(0,int(n/2)+1):
                                    for k in range(0,int((n-2*t)/2)+1):
                                        A3 =  A(n,t,k,az,bz,Ra[2],Rb[2],Rn[2],Rp[2])

                                        nu = l + m + n - 2 * (r + s + t) - (i + j + k)

                                        ff = F(nu,g*np.dot(Rp-Rn,Rp-Rn))

                                        Vn += A1 * A2 * A3 * ff

    Na = norm(ax,ay,az,aa)
    Nb = norm(bx,by,bz,bb)

    Vn *= - Zn * Na * Nb * c * 2 * np.pi / g

    return Vn

def electronic(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,aa,bb,cc,dd,Ra,Rb,Rc,Rd):
    pass
