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

from RHF import *

from matrices import *
from integrals import *
from basis import *
from molecules import *

import numpy as np
import numpy.linalg as la

###########################
###########################
###########################

mol = H2O # Molecule
bs = sto3g_H2O # Basis set
N = 10 # Number of electrons

maxiter = 100 # Maximal number of iteration

verbose = True # Print each SCF step

###########################
###########################
###########################

# Basis set size
K = bs.K

print("Computing overlap matrix S...")
S = S_overlap(bs)

if verbose:
    print(S)

print("Computing orthogonalization matrix X...")
X = X_transform(S)

if verbose:
    print(X)

print("Computing core Hamiltonian...")
Hc = H_core(bs,mol)

if verbose:
    print(Hc)

print("Computing two-electron integrals...")
ee = EE_list(bs)

if verbose:
    print_EE_list(ee)

Pnew = np.zeros((K,K))
P = np.zeros((K,K))

converged = False

print("   ##################")
print("   Starting SCF cycle")
print("   ##################")

iter = 1
while not converged and iter <= maxiter:
    print("\n\n\n#####\nSCF cycle " + str(iter) + ":")
    print("#####")

    Pnew, F, E = RHF_step(bs,mol,N,Hc,X,P,ee,verbose) # Perform an SCF step

    # Print results of the SCF step
    print("\nTotal energy:", energy_tot(P,F,Hc,mol),"\n")
    print("   Orbital energies:")
    print("   ", np.diag(E))

    # Check convergence of the SCF cycle
    if delta_P(P,Pnew) < 1e-12:
        converged = True

        print("\n\n\nTOTAL ENERGY:", energy_tot(P,F,Hc,mol)) # Print final, total energy

    if iter == maxiter:
        print("SCF NOT CONVERGED!")

    P = Pnew

    iter += 1
