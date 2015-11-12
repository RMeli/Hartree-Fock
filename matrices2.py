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
    """

    # H2
    H2 = [Atom("H",(0,0,0),["1s"]),Atom("H",(0,0,1.4),["1s"])]

    # Create the basis set
    sto3g_H2 = STO3G(H2)

    # Overlap matrix
    S_H2 = S_overlap(sto3g_H2)
    T_H2 = T_kinetic(sto3g_H2)

    print("###########")
    print("H2 molecule")
    print("###########")

    print("\nOverlap matrix S:")
    print(S_H2)

    print("\nKinetic matrix T:")
    print(T_H2)

    # HeH+
    HeH = [Atom("H",(0,0,0),["1s"]),Atom("He",(0,0,1.4),["1s"])]

    # H2O
    H2O = [   Atom("H",(0,+1.43233673,-0.96104039),["1s"]),
                Atom("H",(0,-1.43233673,-0.96104039),["1s"]),
                Atom("O",(0,0,0.24026010),["1s","2s","2p"])]

    sto3g_H2O = STO3G(H2O)

    # Overlap matrix
    S_H2O = S_overlap(sto3g_H2O)
    T_H2O = T_kinetic(sto3g_H2O)

    print("\n\n\n")
    print("############")
    print("H2O molecule")
    print("############")

    print("\nOverlap matrix S:")
    print(S_H2O)

    print("\nKinetic matrix T:")
    print(T_H2O)
