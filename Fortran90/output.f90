MODULE OUTPUT

    IMPLICIT NONE

    CONTAINS

        SUBROUTINE append_xyz(Nn,atoms,pos,unit)

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: Nn                               ! Number of particles
            CHARACTER(len=2), dimension(Nn), intent(in) :: atoms    ! Atom type
            REAL*8, dimension(Nn,3), intent(in) :: pos              ! Particles' position
            INTEGER, intent(in) :: unit                             ! Output unit

            ! INTERMEDIATE VARIABLES
            INTEGER :: i                                    ! Loop index

            WRITE(unit,*) Nn
            WRITE(unit,*) ! TODO: add step number
            DO i = 1, Nn
                WRITE(unit,*) atoms(i), pos(i,:)
            END DO

        END SUBROUTINE

END MODULE OUTPUT
