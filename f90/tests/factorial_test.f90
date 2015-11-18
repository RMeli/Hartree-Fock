!-----
! TEST
!-----
PROGRAM test

    USE FACT

    IMPLICIT NONE

    INTEGER :: i
    INTEGER, PARAMETER :: N = 10

    DO i = 0,N
        WRITE(*,*) i, "! =", factorial(i)
    END DO

    DO i = -1,N
        WRITE(*,*) i, "! =", factorial2(i)
    END DO

END PROGRAM test
