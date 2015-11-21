MODULE MOLECULE

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE load_molecule(fname,Nel,Natm,atoms,coordinates)
        IMPLICIT NONE

        INTEGER, intent(out) :: Nel, Natm

        CHARACTER, dimension(Natm,2) :: atoms ! Name of atoms
        REAL*8, dimension(Natm,3), intent(out) :: coordinates ! Coordinates of atoms

    END SUBROUTINE load_molecule

END MODULE MOLECULE
