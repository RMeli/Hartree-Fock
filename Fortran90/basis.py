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


class Atom:
    """
    Class representing an atom.
    """

    def __init__(self,name,R,Z,orbitals):
        """
        Initializer for ATOM

        INPUT:
            NAME: Name of the element
            R: Position (cartesian coordinates, atomic units)
            Z: Atomic charge
            ORBITALS: list of orbitals for this atom
        """

        self.name = name
        self.R = R
        self.orbitals = orbitals
        self.Z = Z

class STO3G():
    """
    STO-3G minimal basis set.
    """
    def __init__(self,atoms):
        """
        Initializer for STO3G

        INPUT:
            ATOMS: list of atoms

        Source

            Modern Quantum Chemistry
            Szabo and Ostlund
            Dover
            1989
        """

        # Exponential coefficients for the Gaussian orbitals
        self.zeta1 = {  "H":1.24,
                        "He":2.0925,
                        "Li":2.69,
                        "Be":3.68,
                        "B":4.68,
                        "C":5.67,
                        "N":6.67,
                        "O":7.66,
                        "F":8.65}
        self.zeta2 = {  "Li":0.75,
                        "Be":1.10,
                        "B":1.45,
                        "C":1.72,
                        "N":1.95,
                        "O":2.25,
                        "F":2.55}

        self.STO3G = []

        for a in atoms: # For every atom
            for o in a.orbitals: # For every atomic orbital
                if o == "1s":
                    a1 = 0.109818 * self.zeta1[a.name]**2
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
                    d3 = 0.1559163

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
        """
        Return the basis set.
        """

        return self.STO3G

    def info(self):
        """
        Print informations about the bais set.
        """
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
