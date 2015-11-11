import numpy as np
import scipy.special as spec

def gaussian_product(aa,bb,Ra,Rb):

    a = aa + bb

    R = (aa * Ra + bb * Rb) / (aa + bb)

    c = np.dot(Ra-Rb,Ra-Rb)
    c *= - aa*bb / (aa + bb)
    c = np.exp(c)

    return a,R,c

def overlap_ss(aa,bb,Ra,Rb):
    a,R,c = gaussian_product(aa,bb,Ra,Rb)

    overlap = np.sqrt((np.pi / a)**3)
    overlap *= c
    overlap *= (2*aa/np.pi)**(3./4.) * (2*bb/np.pi)**(3./4.) # Normalization

    return overlap

def kinetic_ss(aa,bb,Ra,Rb):
    o = overlap_ss(aa,bb,Ra,Rb)

    k = aa * bb / (aa + bb)
    k *= 3 - 2 * k * np.dot(Ra-Rb,Ra-Rb)
    k *= o # Normalization is already contained here

    return k

def nuclear_ss(aa,bb,Ra,Rb,Zn,Rn):
    a,R,c = gaussian_product(aa,bb,Ra,Rb)

    n = -2 * np.pi / a * Zn * c
    n *= F_ss(a * np.dot(R-Rn,R-Rn))
    n *= (2*aa/np.pi)**(3./4.) * (2*bb/np.pi)**(3./4.) # Normalization

    return n

def electronic_ss(aa,bb,cc,dd,Ra,Rb,Rc,Rd):
    p,Rp,cp = gaussian_product(aa,bb,Ra,Rb)
    q,Rq,cq = gaussian_product(cc,dd,Rc,Rd)

    e = 2 * np.pi**(5./2.) / (p * q * np.sqrt(p+q))
    e *= cp * cq
    e *= F_ss(p * q / (p + q) * np.dot(Rp-Rq,Rp-Rq))
    e *= (2*aa/np.pi)**(3./4.) * (2*bb/np.pi)**(3./4.) # Normalization
    e *= (2*cc/np.pi)**(3./4.) * (2*dd/np.pi)**(3./4.) # Normalization

    return e

def EE_list(basis,N,R):
    """
    Multidimensional array of two-electron integrals.
    """

    EE = np.zeros((N,N,N,N))

    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    for a in basis[i]:
                        for b in basis[j]:
                            for c in basis[k]:
                                for d in basis[l]:
                                    f = a[0].conjugate() * b[0].conjugate() * c[0] * d[0]
                                    EE[i,j,k,l] += f * electronic_ss(a[1],b[1],c[1],d[1],R[i],R[j],R[k],R[l])

    return EE

def print_EE_list(ee,N):
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    print("({0},{1},{2},{3})  {4}".format(i+1,j+1,k+1,l+1,ee[i,j,k,l]))

def F_ss(t,cutoff=1e-6):
    if abs(t) < cutoff:
        return 1. - t / 3.
    else:
        # Eq. A.32
        return 0.5 * np.sqrt(np.pi / t) * spec.erf(np.sqrt(t))

if __name__ == "__main__":
    z1 = 2.0925 # He
    z2 = 1.24 # H

    # Basis set
    STO3G = [[(0.444635,0.109818*z1**2),(0.535328,0.405771*z1**2),(0.154329,2.22766*z1**2)],
        [(0.444635,0.109818*z2**2),(0.535328,0.405771*z2**2),(0.154329,2.22766*z2**2)]]

    N = 2 # Number of electrons

    Z = (2,1) # Atomic charges

    R = (0.0,1.4632) # Atomic positions

    print("\nTwo-electron integrals:")
    print_EE_list(EE_list(STO3G,N,R),N)
