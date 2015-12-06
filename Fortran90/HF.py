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

from basis import *

import os
import sys

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

        atom_idx = 0 # Atom index (start from 1)

        for a in atoms: # For every atom

            atom_idx += 1

            for o in a.orbitals: # For every atomic orbital
                if o == "1s":
                    a1 = 0.109818 * self.zeta1[a.name]**2
                    a2 = 0.405771 * self.zeta1[a.name]**2
                    a3 = 2.22766 * self.zeta1[a.name]**2
                    d1 = 0.444635
                    d2 = 0.535328
                    d3 = 0.154329

                    self.STO3G.append({  "AOn":a.name,#
                                    "AOt":o,
                                    "Aidx":atom_idx,#
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
                                "Aidx":atom_idx,
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
                                    "Aidx":atom_idx,
                                    "R":a.R,
                                    "lx":1,
                                    "ly":0,
                                    "lz":0,
                                    "a":(a1,a2,a3),
                                    "d":(d1,d2,d3)})
                    self.STO3G.append({  "AOn":a.name,
                                    "AOt":o,
                                    "Aidx":atom_idx,
                                    "R":a.R,
                                    "lx":0,
                                    "ly":1,
                                    "lz":0,
                                    "a":(a1,a2,a3),
                                    "d":(d1,d2,d3)})
                    self.STO3G.append({  "AOn":a.name,
                                    "AOt":o,
                                    "Aidx":atom_idx,
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

class IO:
    """
    Prepare input file for Fortran program.
    """
    def __init__(self,fname):

        # Open input file
        self.ifile = open(fname,'r')

        # Create name for the output file
        fnameout = os.path.splitext(fname)[0] + "_f.in"

        # Open output file
        self.ofile = open(fnameout, 'w') # Open output file (fortran input file)

        # Occupied orbitals
        self.orbitals = {   "H" : ["1s"],
                            "He" : ["1s"],
                            "C" : ["1s","2s","2p"],
                            "N" : ["1s","2s","2p"],
                            "O" : ["1s","2s","2p"],
                            "F" : ["1s","2s","2p"]}

        # Nuclear charges
        self.charges = {    "H" : 1,
                            "He" : 2,
                            "C" : 6,
                            "N" : 7,
                            "O" : 8,
                            "F" : 9 }

        self.read() # Read from file

    def __del__(self):
        self.ifile.close() # Close file input file
        self.ofile.close() # Close output file

    def read(self):
        self.calculation = self.ifile.readline().strip() # Calculation
        print(self.calculation)

        self.bsname = self.ifile.readline().strip() # Basis set name

        self.atoms = [] # List of atoms composing the system

        self.Ne = int(self.ifile.readline().strip()) # Number of electrons

        if (self.calculation == "RHF" and self.Ne % 2 != 0):
            print("ERROR: Odd number of electrons for a RHF calculation.")
            sys.exit(-1)

        self.Nn = 0 # Number of nuclei

        for line in self.ifile:
            name, rx, ry, rz = line.split()

            a = Atom(name,[rx,ry,rz],self.charges[name],self.orbitals[name])

            self.atoms.append(a)

            self.Nn += 1

        if self.bsname == "STO-3G":
            self.bs = STO3G(self.atoms)
            self.contractions = 3
        else:
            print("ERROR: This basis set is not implemented.")

    def fortran_input(self):
        self.ofile.write(self.calculation + os.linesep) # Calculation
        self.ofile.write(str(self.Ne) + os.linesep) # Number of electrons
        self.ofile.write(str(self.Nn) + os.linesep) # Number of nuclei
        self.ofile.write(str(self.bs.K) + os.linesep) # Number of basis set
        self.ofile.write(str(self.contractions) + os.linesep) # Number of contractions

        for a in self.atoms:
            Rx = "%+.10e" % float(a.R[0])
            Ry = "%+.10e" % float(a.R[1])
            Rz = "%+.10e" % float(a.R[2])

            Z = "%i" % a.Z

            line = Rx + ' ' + Ry + ' ' + Rz + ' ' + Z +os.linesep

            self.ofile.write(line)

        basis = self.bs.basis()

        for b in basis:
            # Basis function center
            Rx = "%+.10e" % float(b["R"][0])
            Ry = "%+.10e" % float(b["R"][1])
            Rz = "%+.10e" % float(b["R"][2])

            # Basis function angular momenta
            lx = "%i" % b["lx"]
            ly = "%i" % b["ly"]
            lz = "%i" % b["lz"]

            line = Rx + ' ' + Ry + ' ' + Rz + ' ' + lx + ' ' + ly + ' ' + lz

            coeff_a = ''

            for ca in b["a"]:
                coeff_a += ' ' + "%+.10e" % float(ca)

            coeff_d = ''

            for cd in b["d"]:
                coeff_d += ' ' + "%+.10e" % float(cd)

            line +=  coeff_a + coeff_d

            line += ' ' + "%i" % b["Aidx"] + os.linesep

            self.ofile.write(line)

        self.ofile.flush() # Ensure the content is written on the file


if __name__ == "__main__":
    fin = sys.argv[1] # Python input file

    io = IO(fin) # Create input/output object
    io.fortran_input() # Create Fortran input file

    ffin = os.path.splitext(fin)[0] + "_f.in" # Fortran input file name

    os.system("./HF.x " + ffin) # Call Fortran Hartree-Fock program

    os.system("rm " + ffin) # Remove input for Fortran program
