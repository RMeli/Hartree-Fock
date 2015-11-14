class Atom:
    def __init__(self,name,R,Z,orbitals):
        self.name = name
        self.R = R
        self.orbitals = orbitals
        self.Z = Z

class STO3G():
    def __init__(self,atoms):

        self.zeta1 = {  "H":1.24,
                        "He":2.0925,
                        "C":5.67,
                        "O":7.66}
        self.zeta2 = {  "C":1.72,
                        "O":2.25}

        self.STO3G = []

        for a in atoms:
            for o in a.orbitals:
                if o == "1s":
                    a1 = 0.109810 * self.zeta1[a.name]**2
                    a2 = 0.405771 * self.zeta1[a.name]**2
                    a3 = 2.22766 * self.zeta1[a.name]**2
                    d1 = 0.444635
                    d2 = 0.535328
                    d3 = 0.154329

                    self.STO3G.append({  "AOn":a.name,#
                                    "AOt":o,#
                                    "R":a.R,#
                                    "lx":0,#
                                    "ly":0,#
                                    "lz":0,#
                                    "a":(a1,a2,a3),#
                                    "d":(d1,d2,d3)})

                if o == "2s":
                    a1 = 0.0751386 * self.zeta2[a.name]**2
                    a2 = 0.231031 * self.zeta2[a.name]**2
                    a3 = 0.994203 * self.zeta2[a.name]**2
                    d1 = 0.700115
                    d2 = 0.399513
                    d3 = -0.0999672

                    self.STO3G.append({  "AOn":a.name,
                                "AOt":o,
                                "R":a.R,
                                "lx":0,
                                "ly":0,
                                "lz":0,
                                "a":(a1,a2,a3),
                                "d":(d1,d2,d3)})

                if o == "2p":
                    a1 = 0.0751386 * self.zeta2[a.name]**2
                    a2 = 0.231031 * self.zeta2[a.name]**2
                    a3 = 0.994203 * self.zeta2[a.name]**2
                    d1 = 0.391957
                    d2 = 0.607684
                    d3 = 0.155916

                    self.STO3G.append({  "AOn":a.name,
                                    "AOt":o,
                                    "R":a.R,
                                    "lx":1,
                                    "ly":0,
                                    "lz":0,
                                    "a":(a1,a2,a3),
                                    "d":(d1,d2,d3)})
                    self.STO3G.append({  "AOn":a.name,
                                    "AOt":o,
                                    "R":a.R,
                                    "lx":0,
                                    "ly":1,
                                    "lz":0,
                                    "a":(a1,a2,a3),
                                    "d":(d1,d2,d3)})
                    self.STO3G.append({  "AOn":a.name,
                                    "AOt":o,
                                    "R":a.R,
                                    "lx":0,
                                    "ly":0,
                                    "lz":1,
                                    "a":(a1,a2,a3),
                                    "d":(d1,d2,d3)})

            self.K = len(self.STO3G)

    def basis(self):
        return self.STO3G

    def info(self):
        print("########################")
        print("STO-3G MINIMAL BASIS SET")
        print("########################\n")

        for b in self.STO3G:
            print(b["AOn"] + " orbital:")
            print("   " + b["AOt"] + ":")
            print("      R = ", b["R"])
            print("      lx = " + str(b["lx"]))
            print("      ly = " + str(b["ly"]))
            print("      lz = " + str(b["lz"]))
            print("      alpha = ", b["a"])
            print("      d = ", b["d"],'\n')


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

    # HeH+
    atoms = [Atom("H",(0,0,0),["1s"]),Atom("He",(0,0,1.4),["1s"])]

    # Create the basis set
    sto3g = STO3G(atoms)

    # Display informations
    sto3g.info()

    print("\n\n\n\n\n")

    # H2O
    atoms = [   Atom("H",(0,+1.43233673,-0.96104039),["1s"]),
                Atom("H",(0,-1.43233673,-0.96104039),["1s"]),
                Atom("O",(0,0,0.24026010),["1s","2s","2p"])]

    # Create the basis set
    sto3g = STO3G(atoms)

    # Display informations
    sto3g.info()
