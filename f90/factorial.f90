MODULE FACT

    IMPLICIT NONE

CONTAINS

    FUNCTION factorial(n) result(fact)
        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: n

        ! INTERMEDIATE VARIABLE
        INTEGER :: i

        ! OUTPUT
        INTEGER :: fact

        fact = 1

        DO i = 2, n
            fact = fact * i
        END DO

    END FUNCTION factorial





    FUNCTION factorial2(n) result(fact2)
        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: n

        ! INTERMEDIATE VARIABLE
        INTEGER :: i

        ! OUTPUT
        INTEGER :: fact2

        IF (n == -1) THEN
            fact2 = -1
        ELSE IF (MOD(n,2) == 0) THEN ! EVEN NUMBER
            fact2 = 1

            DO i = 2, n, 2
                fact2 = fact2 * i
            END DO
        ELSE ! ODD NUMBER
            fact2 = 1

            DO i = 1, n, 2
                fact2 = fact2 * i
            END DO
        END IF

    END FUNCTION factorial2

END MODULE FACT




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
