from basis2 import *
from integrals2 import *

def S_overlap(basis):
    """
    Compute overlap matrix S.

    BASIS: basis set.
    """
    # Size of the basis set
    K = basis.K

    # List of basis functions
    B = basis.basis()

    S = np.zeros((K,K))

    for i,b1 in enumerate(B):
        for j,b2 in enumerate(B):
            for a1,d1 in zip(b1["a"],b1["d"]):
                for a2,d2 in zip(b2["a"],b2["d"]):
                    R1 = b1["R"]
                    R2 = b2["R"]

                    tmp = d1.conjugate()*d2
                    tmp *= overlap(b1["lx"],b1["ly"],b1["lz"],b2["lx"],b2["ly"],b2["lz"],a1,a2,R1,R2)

                    S[i,j] +=  tmp

    return S

def T_kinetic(basis):
    """
    Compute kinetic matrix T.

    BASIS: basis set.
    """
    # Size of the basis set
    K = basis.K

    # List of basis functions
    B = basis.basis()

    T = np.zeros((K,K))

    for i,b1 in enumerate(B):
        for j,b2 in enumerate(B):
            for a1,d1 in zip(b1["a"],b1["d"]):
                for a2,d2 in zip(b2["a"],b2["d"]):
                    R1 = b1["R"]
                    R2 = b2["R"]

                    tmp = d1.conjugate()*d2
                    tmp *= kinetic(b1["lx"],b1["ly"],b1["lz"],b2["lx"],b2["ly"],b2["lz"],a1,a2,R1,R2)

                    T[i,j] +=  tmp

    return T

def V_nuclear(basis,atom):
    """
    Compute nuclear-electron potential energy matrix Vn.

    BASIS: basis set.
    RN: Nuclear position.
    ZN: Nuclear charge.
    """
    # Size of the basis set
    K = basis.K

    # List of basis functions
    B = basis.basis()

    # Nuclear coordinates
    Rn = atom.R

    # Nuclear charge
    Zn = atom.Z


    Vn = np.zeros((K,K))

    for i,b1 in enumerate(B):
        for j,b2 in enumerate(B):
            for a1,d1 in zip(b1["a"],b1["d"]):
                for a2,d2 in zip(b2["a"],b2["d"]):
                    R1 = b1["R"]
                    R2 = b2["R"]

                    tmp = d1.conjugate()*d2
                    tmp *= nuclear(b1["lx"],b1["ly"],b1["lz"],b2["lx"],b2["ly"],b2["lz"],a1,a2,R1,R2,Rn,Zn)

                    Vn[i,j] +=  tmp

    return Vn

def H_core(basis,molecule):
    T = T_kinetic(basis)

    # Size of the basis set
    K = basis.K

    Vn = np.zeros((K,K))

    for atom in molecule:
        Vn += V_nuclear(basis,atom)

    return T + Vn

if __name__ == "__main__":

    """
    Results compared with

        Modern Quantum Chemistry
        Szabo and Ostlund
        Dover
        1989

    and

        The Mathematica Journal
        Evaluation of Gaussian Molecular Integrals
        I. Overlap Integrals
        Minhhuy Hô and Julio Manuel Hernández-Pérez
        2012

    and

        The Mathematica Journal
        Evaluation of Gaussian Molecular Integrals
        II. Kinetic-Energy Integrals
        Minhhuy Hô and Julio Manuel Hernández-Pérez
        2013

    and

        The Mathematica Journal
        Evaluation of Gaussian Molecular Integrals
        III. Nuclear-Electron attraction Integrals
        Minhhuy Hô and Julio Manuel Hernández-Pérez
        2014
    """

    # H2
    H2 = [Atom("H",(0,0,0),1,["1s"]),Atom("H",(0,0,1.4),1,["1s"])]

    # Create the basis set
    sto3g_H2 = STO3G(H2)

    # Compute matrices
    S_H2 = S_overlap(sto3g_H2)
    T_H2 = T_kinetic(sto3g_H2)
    Vn1_H2 = V_nuclear(sto3g_H2,H2[0])
    Vn2_H2 = V_nuclear(sto3g_H2,H2[1])
    H_core_H2 = H_core(sto3g_H2,H2)

    print("###########")
    print("H2 molecule")
    print("###########")

    print("\nOverlap matrix S:")
    print(S_H2)

    print("\nKinetic matrix T:")
    print(T_H2)

    print("\nElectron-nucleus interaction " + H2[0].name + " :")
    print(Vn1_H2)

    print("\nElectron-nucleus interaction " + H2[1].name + " :")
    print(Vn2_H2)

    print("\nCore Hamiltonian:")
    print(H_core_H2)

    # HeH+
    HeH = [Atom("H",(0,0,0),1,["1s"]),Atom("He",(0,0,1.4632),2,["1s"])]

    # Create the basis set
    sto3g_HeH = STO3G(HeH)

    # Compute matrices
    S_HeH = S_overlap(sto3g_HeH)
    T_HeH = T_kinetic(sto3g_HeH)
    Vn1_HeH = V_nuclear(sto3g_HeH,HeH[0])
    Vn2_HeH = V_nuclear(sto3g_HeH,HeH[1])
    H_core_HeH = H_core(sto3g_HeH,HeH)

    print("\n\n\n")
    print("############")
    print("HeH molecule")
    print("############")

    print("\nOverlap matrix S:")
    print(S_HeH)

    print("\nKinetic matrix T:")
    print(T_HeH)

    print("\nElectron-nucleus interaction " + HeH[0].name + " :")
    print(Vn1_HeH)

    print("\nElectron-nucleus interaction " + HeH[1].name + " :")
    print(Vn2_HeH)

    print("\nCore Hamiltonian:")
    print(H_core_HeH)

    # H2O
    H2O = [   Atom("H",(0,+1.43233673,-0.96104039),1,["1s"]),
                Atom("H",(0,-1.43233673,-0.96104039),1,["1s"]),
                Atom("O",(0,0,0.24026010),8,["1s","2s","2p"])]

    sto3g_H2O = STO3G(H2O)

    # Overlap matrix
    S_H2O = S_overlap(sto3g_H2O)
    T_H2O = T_kinetic(sto3g_H2O)
    Vn1_H2O = V_nuclear(sto3g_H2O,H2O[0])
    Vn2_H2O = V_nuclear(sto3g_H2O,H2O[1])
    Vn3_H2O = V_nuclear(sto3g_H2O,H2O[2])

    Vn_H2O = Vn1_H2O + Vn2_H2O + Vn3_H2O

    H_core_H2O = H_core(sto3g_H2O,H2O)

    print("\n\n\n")
    print("############")
    print("H2O molecule")
    print("############")

    print("\nOverlap matrix S:")
    print(S_H2O)

    print("\nKinetic matrix T:")
    print(T_H2O)

    print("\nTotal electron-nucleus interaction:")
    print(Vn_H2O)

    print("\nCore hamiltonian:")
    print(H_core_H2O)
