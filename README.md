# Hartree-Fock 
 
Hartree-Fock code for total energy calculations within Pople's STO-3G basis set. 
 
Other STO-kG basis sets can be easily added. The general implementation follows reference [1].  
One- and two-electron integrals were originally implemented using analytical results as explained in Reference [2] (very slow). Overlap and multiplole integrals are now implemented using the Obara-Saika scheme as presented in Regerence [3]. 
 
## Python 
 
First iteration of the program. Restructed Hartree-Fock (RHF) calculations only. Extremely slow. 
 
## Fortran90 
 
Second iteration of the program. RHF and Unrestricted Hartree-Fock (UHF). Much faster, but only because of compiler optimization. 
 
Total energy derivatives with respect to nuclear coordinates are computed using finite differences (very unefficient). Born-Oppenheimer molecular dynamics (BOMD) is not working. 
 
### Usage 
 
 
## Known Bugs 
 
- BOMD is broken. 
 
## TODO 
 
- Reduce libraries requirements. 
- Faster implementation of one- and two-electron integrals. 
- Analytical derivatives. 
- Fix BOMD. 
- Visualization of molecular orbitals. 
 
## References 
 
[1] A. Szabo and N. S. Ostlund, *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory*, Dover Publications (1996). 
 
[2] D. Cook, *Handbook of Computational Chemistry*, Oxford University Press (1998). 
 
[3] T. Helgaker, P. Jørgensen and J. Olsen, *Molecular Electronic-Structure Theory*, Wiley (2000).*
