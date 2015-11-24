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

mol = HeH # Molecule
bs = sto3g_HeH # Basis set
N = 2 # Number of electrons

maxiter = 2

verbose = True

###########################
###########################
###########################

# Basis set size
K = bs.K

print("Computing overlap matrix S...")
S = S_overlap(bs)

print("Computing orthogonalization matrix X...")
X = X_transform(S)

print("Computing core Hamiltonian...")
Hc = H_core(bs,mol)

print("Computing two-electron integrals...")
ee = EE_list(bs)

print_EE_list(ee)

Pnew = np.zeros((K,K))
P = np.zeros((K,K))

converged = False

print("   ##################")
print("\n\n\nStarting SCF cycle")
print("   ##################")

iter = 1
while not converged and iter <= maxiter:
    print("\n\n\n#####\nSCF cycle " + str(iter) + ":")
    print("#####")

    Pnew, F, E = RHF_step(bs,mol,N,Hc,X,P,ee,verbose)

    print("\nTotal energy:", energy_tot(P,F,Hc,mol),"\n")
    print("   Orbital energies:")
    print("   ", np.diag(E))

    if delta_P(P,Pnew) < 1e-12:
        converged = True

        print("\n\n\nTOTAL ENERGY:", energy_tot(P,F,Hc,mol))

    if iter == maxiter:
        print("SCF NOT CONVERGED!")

    P = Pnew

    iter += 1
