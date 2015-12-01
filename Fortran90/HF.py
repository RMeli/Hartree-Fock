from basis import *

import os

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
                            "N" : ["1s","2s","2p"],
                            "O" : ["1s","2s","2p"]}

        # Nuclear charges
        self.charges = {    "H" : 1,
                            "He" : 2,
                            "N" : 7,
                            "O" : 8 }

        self.read() # Read from file

    def __del__(self):
        pass
        #self.ifile.close() # Close file input file
        #self.ofile.close() # Close output file

    def read(self):
        self.bsname = self.ifile.readline().strip() # Basis set name

        self.atoms = [] # List of atoms composing the system

        self.Ne = self.ifile.readline().strip() # Number of electrons

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
        self.ofile.write(str(self.Ne) + os.linesep) # Number of electrons
        self.ofile.write(str(self.Nn) + os.linesep) # Number of nuclei
        self.ofile.write(str(self.bs.K) + os.linesep) # Number of basis set
        self.ofile.write(str(self.contractions) + os.linesep) # Number of contractions

        for a in self.atoms:
            Rx = "%+.10e" % float(a.R[0])
            Ry = "%+.10e" % float(a.R[1])
            Rz = "%+.10e" % float(a.R[2])

            line = Rx + ' ' + Ry + ' ' + Rz + os.linesep

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

            line +=  coeff_a + coeff_d + os.linesep


            self.ofile.write(line)

        self.ofile.flush() # Ensure the content is written on the file


if __name__ == "__main__":
    import sys

    io = IO(sys.argv[1]) # Open file as argument
    io.fortran_input()
