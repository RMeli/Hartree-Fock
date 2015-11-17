from matrices import *

import numpy.linalg as la

def RHF_step(basis,molecule,N,H,X,P_old,ee,verbose=False):

    if verbose:
        print("\nDensity matrix P:")
        print(P_old)

    G = G_ee(basis,molecule,P_old,ee)

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

    # Sort eigenvalues from smallest to highest (needed to compute P correctly)
    idx = e.argsort()
    e = e[idx]
    Cx = Cx[:,idx]

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
        print(C)

    Pnew = P_density(C,N)

    return Pnew, F, e

def delta_P(P_old,P_new):
    delta = 0

    n = P_old.shape[0]

    for i in range(n):
        for j in range(n):
            delta += (P_old[i,j] - P_new[i,j])**2

    return (delta / 4.)**(0.5)

def energy_el(P,F,H):

    # Size of the basis set
    K = P.shape[0]

    E = 0

    for i in range(K):
        for j in range(K):
            E += 0.5 * P[i,j] * (H[i,j] + F[i,j])

    return E

def energy_n(molecule):

    en = 0

    for i in range(len(molecule)):
        for j in range(i+1,len(molecule)):
            atomi = molecule[i]
            atomj = molecule[j]

            Ri = np.asarray(atomi.R)
            Rj = np.asarray(atomj.R)

            en += atomi.Z * atomj.Z / la.norm(Ri - Rj)

    return en

def energy_tot(P,F,H,molecule):
    return energy_el(P,F,H) + energy_n(molecule)
