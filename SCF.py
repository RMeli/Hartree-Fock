from matrices import *

import numpy.linalg as la

def SCF_RHF_step(basis,N,R,Z,H,X,P_old,ee,verbose=False):

    if verbose:
        print("\nDensity matrix P:")
        print(P_old)

    G = G_ee(basis,N,R,Z,P_old,ee)

    if verbose:
        print("\nG matrix:")
        print(G)

    F = H + G

    if verbose:
        print("\nFock matrix:")
        print(F)

    Fx = np.dot(X.conj().T,np.dot(F,X))

    if verbose:
        print("\nFock matrix in orthogonal orbital basis:")
        print(Fx)

    e, Cx = la.eigh(Fx)

    if verbose:
        print("\nCoefficients in orthogonal orbital basis:")
        print(Cx)

    e = np.diag(e)

    if verbose:
        print("\nEnergies in orthogonal orbital basis:")
        print(e)

    C = np.dot(X,Cx)

    if verbose:
        print("\nCoefficients:")
        print(e)

    Pnew = P_density(C,N)

    return Pnew, F, H

def delta_P(P_old,P_new):
    delta = 0

    n = P_old.shape[0]

    for i in range(n):
        for j in range(n):
            delta += (P_old[i,j] - P_new[i,j])**2

    return (delta / 4.)**(0.5)

def energy_el(P,F,H,N):
    E = 0

    for i in range(N):
        for j in range(N):
            E += 0.5 * P[i,j] * (H[i,j] + F[i,j])

    return E

def energy_n(Z,R):
    en = 0

    for i in range(len(R)):
        for j in range(i+1,len(R)):
            en += - Z[i] * Z[j] / (R[i] - R[j])

    return en

def energy_tot(P,F,H,N,Z,R):
    return energy_el(P,F,H,N) + energy_n(Z,R)
