from SCF import *
from integrals import *
from basis import *

import numpy.linalg as la

# Basis set
z1 = 2.09 # He
z2 = 1.24 # H



N = 2 # Number of electrons
Z = (1,1) # Atomic chanrges
R = (0.0,1.4) # Atomic distances

STO3G = STO3G_basis([1.24]*N)

S = S_overlap(STO3G,N,R) # Overlap matrix
H = H_core(STO3G,N,R,Z) # Core Hamiltonian
ee = EE_list(STO3G,N,R) # Two-electron integrals
X = X_transform(S)


Pnew = np.zeros((N,N))
P = np.zeros((N,N))

converged = False

i = 1
while not converged:
    print("SCF cycle " + str(i) + ":")

    Pnew, F, H = SCF_RHF_step(STO3G,N,R,Z,H,X,P,ee)

    print("Total energy:", energy_tot(P,F,H,N,Z,R),"\n\n\n")

    if delta_P(P,Pnew) < 1e-12:
        converged = True

        print("TOTAL ENERGY:", energy_tot(P,F,H,N,Z,R),"\n\n\n")

    P = Pnew

    i += 1
