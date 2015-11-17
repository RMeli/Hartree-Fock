from RHF import *

from matrices2 import *
from integrals2 import *
from basis2 import *

import numpy as np
import numpy.linalg as la

# He
He = [Atom("He",(0,0,0),2,["1s"])]
sto3g_He = STO3G(He)

# B
Be = [Atom("Be",(0,0,0),4,["1s","2s"])]
sto3g_Be = STO3G(Be)

# N2
N2 = [  Atom("N",(0,0,0),7,["1s","2s","2p"]),
        Atom("N",(0,0,2.074),7,["1s","2s","2p"])]
sto3g_N2 = STO3G(N2)

# HF
HF = [  Atom("H",(0,0,0),1,["1s"]),
        Atom("F",(0,0,1.807),9,["1s","2s","2p"])]
sto3g_HF = STO3G(HF)

# B
O = [Atom("O",(0,0,0),8,["1s","2s","2p"])]
sto3g_O = STO3G(O)

# H2
H2 = [Atom("H",(0,0,0),1,["1s"]),Atom("H",(0,0,1.4),1,["1s"])]
sto3g_H2 = STO3G(H2)

# HeH+
HeH = [Atom("He",(0,0,1.4632),2,["1s"]),Atom("H",(0,0,0),1,["1s"])]
sto3g_HeH = STO3G(HeH)

# H2O
H2O = [ Atom("H",(1.638036840407,1.136548822547,0.000000000000),1,["1s"]),
        Atom("H",(-1.638036840407,1.136548822547,0.000000000000),1,["1s"]),
        Atom("O",(0.000000000000,-0.143225816552,0.000000000000),8,["1s","2s","2p"])]
sto3g_H2O = STO3G(H2O)


# Molecule
mol = H2O

# Basis set
bs = sto3g_H2O

# Number of electrons
N = 10

# Basis set size
K = bs.K

print("Computing overlap matrix S...")
S = S_overlap(bs)

print(S)

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

maxiter = 50

iter = 1
while not converged and iter <= maxiter:
    print("\n\n\n#####\nSCF cycle " + str(iter) + ":")
    print("#####")

    Pnew, F, E = RHF_step(bs,mol,N,Hc,X,P,ee,True)

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
