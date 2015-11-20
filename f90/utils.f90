MODULE UTILS

    IMPLICIT NONE

    CONTAINS

    ! ------------
    ! PRINT MATRIX
    ! ------------
    SUBROUTINE print_real_matrix(r,c,M)

        ! INPUT
        INTEGER, intent(in) :: r, c
        REAL*8, dimension(r,c), intent(in) :: M

        ! VARIABLES
        INTEGER :: i

        DO i = 1, r
            WRITE(*,'(20G12.6)') M(i,1:c)
        END DO

    END SUBROUTINE print_real_matrix

END MODULE UTILS
