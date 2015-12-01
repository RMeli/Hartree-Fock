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

from matrices import *

import numpy.linalg as la

def RHF_step(basis,molecule,N,H,X,P_old,ee,verbose=False):
    """
    Restricted Hartree-Fock self-consistent field setp.

    INPUT:
        BASIS: basis set
        MOLECULE: molecule, collection of atom object
        N: Number of electrons
        H: core Hamiltonian
        X: tranformation matrix
        P_OLD: Old density matrix
        EE: List of electron-electron Integrals
        VERBOSE: verbose flag (set True to print everything on screen)
    """

    if verbose:
        print("\nDensity matrix P:")
        print(P_old)

    G = G_ee(basis,molecule,P_old,ee) # Compute electron-electron interaction matrix

    if verbose:
        print("\nG matrix:")
        print(G)

    F = H + G # Compute Fock matrix

    if verbose:
        print("\nFock matrix:")
        print(F)

    Fx = np.dot(X.conj().T,np.dot(F,X)) # Compute Fock matrix in the orthonormal basis set (S=I in this set)

    if verbose:
        print("\nFock matrix in orthogonal orbital basis:")
        print(Fx)

    e, Cx = la.eigh(Fx) # Compute eigenvalues and eigenvectors of the Fock matrix

    # Sort eigenvalues from smallest to highest (needed to compute P correctly)
    idx = e.argsort()
    e = e[idx]
    Cx = Cx[:,idx]

    if verbose:
        print("\nCoefficients in orthogonal orbital basis:")
        print(Cx)

    e = np.diag(e) # Extract orbital energies as vector

    if verbose:
        print("\nEnergies in orthogonal orbital basis:")
        print(e)

    C = np.dot(X,Cx) # Transform coefficient matrix in the orthonormal basis to the original basis

    if verbose:
        print("\nCoefficients:")
        print(C)

    Pnew = P_density(C,N) # Compute the new density matrix

    return Pnew, F, e

def delta_P(P_old,P_new):
    """
    Compute the difference between two density matrices.

    INTPUT:
        P_OLD: Olde density matrix
        P_NEW: New density matrix
    OUTPUT:
        DELTA: difference between the two density matrices

    Source:
        Modern Quantum Chemistry
        Szabo and Ostlund
        Dover
        1989
    """
    delta = 0

    n = P_old.shape[0]

    for i in range(n):
        for j in range(n):
            delta += (P_old[i,j] - P_new[i,j])**2

    return (delta / 4.)**(0.5)

def energy_el(P,F,H):
    """
    Compute electronic energy.

    INPUT:
        P: Density matrix
        F: Fock matrix
        H: Core Hamiltonian

    Source:
        Modern Quantum Chemistry
        Szabo and Ostlund
        Dover
        1989
    """

    # Size of the basis set
    K = P.shape[0]

    E = 0

    for i in range(K):
        for j in range(K):
            E += 0.5 * P[i,j] * (H[i,j] + F[i,j])

    return E

def energy_n(molecule):
    """
    Compute nuclear energy (classical nucleus-nucleus repulsion)

    INPUT:
        MOLECULE: molecule, as a collection of atoms
    OUTPUT:
        ENERGY_N: Nuclear energy
    """

    en = 0

    for i in range(len(molecule)):
        for j in range(i+1,len(molecule)):
            # Select atoms from molecule
            atomi = molecule[i]
            atomj = molecule[j]

            # Extract distance from atom
            Ri = np.asarray(atomi.R)
            Rj = np.asarray(atomj.R)

            en += atomi.Z * atomj.Z / la.norm(Ri - Rj)

    return en

def energy_tot(P,F,H,molecule):
    """
    Compute total energy (electronic plus nuclear).

    INPUT:
        P: Density matrix
        F: Fock matrix
        H: Core Hamiltonian
        MOLECULE: molecule, as a collection of atoms
    OUTPUT:
        ENERGY_TOT: total energy
    """
    return energy_el(P,F,H) + energy_n(molecule)
