# Hartree-Fock 
 
Hartree-Fock code for total energy calculations within Pople's STO-3G basis set. 
 

### Details

The general implementation follows Reference [1]. One- and two-electron integrals were originally implemented using analytical results as explained in Reference [2] (very slow). Overlap and multipole integrals are now implemented using the Obara-Saika scheme as presented in Reference [3]. 

### Disclaimer

This program was written in 2015 for the class *Computational Methods in Molecular Quantum Mechanics* at *Ecole Polytechnique Fédérale de Lausanne* and was awarded the maximum grade. 

It was my very first Fortran program, therefor the coding style might be completely off. My Python skill were also very much in development back then.
 
## Python 
 
First iteration of the program. Restricted Hartree-Fock (RHF) calculations only. Extremely slow. 
 
## Fortran90 
 
Second iteration of the program. RHF and Unrestricted Hartree-Fock (UHF). Much faster, because of compiler optimization. 
 
Total energy derivatives with respect to nuclear coordinates are computed only using finite differences (no analytical derivatives, very inefficient). Born-Oppenheimer molecular dynamics (BOMD) is not working properly. 
 
 
## Known Bugs 
 
- BOMD is not working
 
## TODO 
 
- Reduce libraries requirements
- Faster implementation of one- and two-electron integrals 
- Analytical derivatives 
- Fix BOMD
- Visualization of molecular orbitals
 
## References 
 
[1] A. Szabo and N. S. Ostlund, *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory*, Dover Publications (1996). 
 
[2] D. Cook, *Handbook of Computational Chemistry*, Oxford University Press (1998). 
 
[3] T. Helgaker, P. Jørgensen and J. Olsen, *Molecular Electronic-Structure Theory*, Wiley (2000).*
